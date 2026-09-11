import argparse
import gzip
from collections import namedtuple
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from itertools import zip_longest
from pathlib import Path
import time

from fuzzysearch import find_near_matches
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser(
        description="Preprocess FASTQ files", add_help=False
    )
    parser.add_argument("--help", action="help", help="Show this help message and exit")
    parser.add_argument("--reads1", required=True, help="Path to R1 FASTQ")
    parser.add_argument("--reads2", required=True, help="Path to R2 FASTQ")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--core", type=int, default=8, help="Number of worker processes")
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument("--compression_level", type=int, default=6, help="Gzip compression level")
    parser.add_argument("--batch_size", type=int, default=DEFAULT_BATCH_SIZE, help="Read pairs per worker batch")
    parser.add_argument("--gzip_output", type=str_to_bool, default=True, help="Compress extracted FASTQs")
    parser.add_argument("--linker1", default="GTGGCCGATGTTTCGCATCGGCGTACGACT", help="Linker 1 sequence")
    parser.add_argument("--linker2", default="ATCCACGTGCTTGAGAGGCCAGAGCATTCG", help="Linker 2 sequence")
    parser.add_argument("--mm_rate", type=float, default=0.05, help="Linker mismatch rate")
    parser.add_argument("--barcodeA_whitelist", help="Barcode A whitelist")
    parser.add_argument("--barcodeB_whitelist", help="Barcode B whitelist")
    parser.add_argument("--correct_barcode", type=str_to_bool, default=False, help="Correct barcodes")
    parser.add_argument("--cb_len", type=int, required=True, help="Concatenated cell-barcode length")
    parser.add_argument("--umi_len", type=int, required=True, help="UMI length")
    return parser.parse_args()


BARCODE_LEN = None
UMI_LEN = None
DEFAULT_BATCH_SIZE = 50000
Match = namedtuple("Match", ["start", "end"])
BatchResult = namedtuple("BatchResult", ["r1_records", "r2_records", "n_reads", "n_reads_passed", "exact_match_stats", "fuzzy_match_stats"])

WORKER_CONTEXT = None


def open_text(path, mode="rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode, encoding="utf-8")


def read_nonempty_lines(path) -> list[str]:
    with open_text(path) as handle:
        return [line.strip() for line in handle if line.strip()]


def iter_fastq(handle):
    while True:
        id_line = handle.readline()
        if not id_line:
            break
        sequence_line = handle.readline()
        plus_line = handle.readline()
        quality_line = handle.readline()
        if not (sequence_line and plus_line and quality_line):
            raise ValueError("Incomplete FASTQ record encountered.")
        if not id_line.startswith("@") or not plus_line.startswith("+"):
            raise ValueError("Invalid FASTQ structure (missing @ or + line).")
        read_id = id_line[1:].strip()
        sequence = sequence_line.strip()
        quality = quality_line.strip()
        if len(sequence) != len(quality):
            raise ValueError(
                f"Length mismatch (seq {len(sequence)} vs qual {len(quality)}) "
                f"at read {read_id}"
            )
        yield read_id, sequence, quality


def str_to_bool(value):
    """Convert str to bool"""
    if isinstance(value, bool):
        return value
    if value.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif value.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError(f'Boolean value expected, got: {value}')


def get_mm_dist(seq, rate, n = 2):
    return max(int(len(seq) * rate), n)


def read_whitelist(whitelist_path):
    return set(read_nonempty_lines(whitelist_path))


def build_barcode_correction_map(whitelist):
    """Map exact and unambiguous single-substitution barcodes to the whitelist."""
    correction_map = {}
    for bc in whitelist:
        correction_map[bc] = bc
    for correct in whitelist:
        for pos in range(len(correct)):
            orig = correct[pos]
            for base in 'ATCG':
                if base == orig:
                    continue
                err = correct[:pos] + base + correct[pos+1:]
                if err not in correction_map:
                    correction_map[err] = correct
                elif correction_map[err] != correct:
                    correction_map[err] = None
    # drop ambiguous
    return {k: v for k, v in correction_map.items() if v is not None}


def iter_paired_fastq_batches(r1_handle, r2_handle, batch_size):
    batch = []
    for r1_record, r2_record in zip_longest(iter_fastq(r1_handle), iter_fastq(r2_handle)):
        if r1_record is None or r2_record is None:
            raise ValueError("R1 and R2 FASTQ files contain different numbers of records.")
        batch.append((r1_record, r2_record))
        if len(batch) >= batch_size:
            yield batch
            batch = []
    if batch:
        yield batch


def find_exact_matches(seq_str, pattern):
    """Return one exact match in the full sequence, or [] if none/ambiguous."""
    first = seq_str.find(pattern)
    if first == -1:
        return []
    second = seq_str.find(pattern, first + 1)
    if second != -1:
        return []
    return [Match(first, first + len(pattern))]


def find_fuzzy_matches(seq_str, pattern, max_errors):
    """Return all fuzzy matches in the full sequence."""
    matches = find_near_matches(pattern, seq_str, max_l_dist=max_errors)
    return [Match(match.start, match.end) for match in matches]


class MatchResult:
    """Container for match results"""
    def __init__(self):
        self.linker1_matches = []
        self.linker2_matches = []
        self.match_method = ""  # "exact", "fuzzy", "mixed", or "failed"
        self.match_stats = [-1, -1, -1]  # [all, linker1, linker2], 0 represents fuzzy, 1 represents exact


def find_all_matches(seq_str, linker1, linker2, linker1_mm, linker2_mm):
    """
    Try exact matches first; fall back to fuzzy per element if needed.
    Search the full read for each linker.
    Require exactly one hit for each element.
    """
    result = MatchResult()

    mixed_success = True
    methods_used = []

    result.linker2_matches = find_exact_matches(seq_str, linker2)
    if len(result.linker2_matches) == 1:
        methods_used.append('exact')
        result.match_stats[2] = 1
    else:
        result.linker2_matches = find_fuzzy_matches(seq_str, linker2, linker2_mm)
        methods_used.append('fuzzy')
        result.match_stats[2] = 0
        if len(result.linker2_matches) != 1:
            mixed_success = False

    if mixed_success:
        result.linker1_matches = find_exact_matches(seq_str, linker1)
        if len(result.linker1_matches) == 1:
            methods_used.append('exact')
            result.match_stats[1] = 1
        else:
            result.linker1_matches = find_fuzzy_matches(seq_str, linker1, linker1_mm)
            methods_used.append('fuzzy')
            result.match_stats[1] = 0
            if len(result.linker1_matches) != 1:
                mixed_success = False

    if mixed_success:
        result.match_method = "mixed" if 'fuzzy' in methods_used else "exact"
        if result.match_stats[1] == 1 and result.match_stats[2] == 1:
            result.match_stats = [1, -1, -1]
        elif result.match_stats[1] == 0 and result.match_stats[2] == 0:
            result.match_stats = [0, -1, -1]
    else:
        result.match_method = "failed"
        result.linker1_matches = []
        result.linker2_matches = []
        result.match_stats = [-1, -1, -1]

    return result


def extract_barcode(seq, qual, linker2_start, linker2_end, linker1_end):
    """
    seq/qual are strings. Extract:
      - barcodeA: configured component length immediately after linker2
      - barcodeB: configured component length immediately before linker2
      - umi: configured length immediately after linker1
    """
    barcodeA = seq[linker2_end: linker2_end + BARCODE_LEN]
    barcodeB = seq[linker2_start - BARCODE_LEN:linker2_start]
    barcodeA_q = qual[linker2_end: linker2_end + BARCODE_LEN]
    barcodeB_q = qual[linker2_start - BARCODE_LEN:linker2_start]
    barcode = barcodeB + barcodeA
    barcode_q = barcodeB_q + barcodeA_q
    umi = seq[linker1_end: linker1_end + UMI_LEN]
    umi_q = qual[linker1_end: linker1_end + UMI_LEN]
    return barcode, umi, barcode_q, umi_q


def correct_barcode_with_correction_map(barcode, barcodeA_correction_map, barcodeB_correction_map):
    """Use global maps built from whitelists; return corrected barcode or None."""
    bcB = barcode[:BARCODE_LEN]
    bcA = barcode[BARCODE_LEN:]
    bcA_cor = barcodeA_correction_map.get(bcA, None)
    if bcA_cor is None:
        return None
    bcB_cor = barcodeB_correction_map.get(bcB, None)
    if bcB_cor is None:
        return None
    return bcB_cor + bcA_cor


def qual_to_string(qual):
    """Compatibility shim: keep as-is when string; (still supports list[int] if ever passed)."""
    if isinstance(qual, str):
        return qual
    return ''.join(chr(q + 33) for q in qual)


def format_seqrecord_fastq(record_id, seq, qual):
    qual_str = qual_to_string(qual)
    return f"@{record_id}\n{seq}\n+\n{qual_str}\n"


class MatchConfig:
    # extract config
    def __init__(self, linker1, linker2, mm_rate):
        self.linker1 = linker1
        self.linker2 = linker2
        self.mm_rate = mm_rate


class BarcodeConfig:
    # barcode correction config
    def __init__(self, barcodeA_whitelist, barcodeB_whitelist):
        self.barcodeA_whitelist = barcodeA_whitelist
        self.barcodeB_whitelist = barcodeB_whitelist


def init_worker(match_config, linker1_mm, linker2_mm, correct_barcode, barcodeA_correction_map, barcodeB_correction_map, barcode_len, umi_len):
    global WORKER_CONTEXT, BARCODE_LEN, UMI_LEN
    BARCODE_LEN = barcode_len
    UMI_LEN = umi_len
    WORKER_CONTEXT = {
        "match_config": match_config,
        "linker1_mm": linker1_mm,
        "linker2_mm": linker2_mm,
        "correct_barcode": correct_barcode,
        "barcodeA_correction_map": barcodeA_correction_map,
        "barcodeB_correction_map": barcodeB_correction_map,
    }


def process_read_batch(batch, context=None):
    if context is None:
        context = WORKER_CONTEXT
    match_config = context["match_config"]
    linker1_mm = context["linker1_mm"]
    linker2_mm = context["linker2_mm"]
    correct_barcode = context["correct_barcode"]
    barcodeA_correction_map = context["barcodeA_correction_map"]
    barcodeB_correction_map = context["barcodeB_correction_map"]

    exact_match_stats = [0, 0, 0]
    fuzzy_match_stats = [0, 0, 0]
    n_reads = 0
    n_reads_passed = 0
    r1_records = []
    r2_records = []

    for (r1_id, r1_seq, r1_qual), (r2_id, r2_seq, r2_qual) in batch:
        n_reads += 1
        match_result = find_all_matches(r1_seq, match_config.linker1, match_config.linker2, linker1_mm, linker2_mm)
        if (len(match_result.linker2_matches) == 1) and (len(match_result.linker1_matches) == 1):
            linker2_start = match_result.linker2_matches[0].start
            linker2_end = match_result.linker2_matches[0].end
            linker1_end = match_result.linker1_matches[0].end
            barcode, umi, barcode_q, umi_q = extract_barcode(r1_seq, r1_qual, linker2_start, linker2_end, linker1_end)
            if correct_barcode:
                barcode = correct_barcode_with_correction_map(barcode, barcodeA_correction_map, barcodeB_correction_map)
            if barcode is not None and len(barcode) == 2 * BARCODE_LEN and len(umi) == UMI_LEN:
                n_reads_passed += 1
                r1_records.append(format_seqrecord_fastq(r1_id, barcode + umi, barcode_q + umi_q))
                r2_records.append(format_seqrecord_fastq(r2_id, r2_seq, r2_qual))
        for i in range(3):
            if match_result.match_stats[i] == 0:
                fuzzy_match_stats[i] += 1
            elif match_result.match_stats[i] == 1:
                exact_match_stats[i] += 1

    return BatchResult(r1_records, r2_records, n_reads, n_reads_passed, exact_match_stats, fuzzy_match_stats)


def add_batch_stats(result, exact_match_stats, fuzzy_match_stats):
    for i in range(3):
        exact_match_stats[i] += result.exact_match_stats[i]
        fuzzy_match_stats[i] += result.fuzzy_match_stats[i]


def write_batch_result(result, out_r1, out_r2):
    out_r1.writelines(result.r1_records)
    out_r2.writelines(result.r2_records)


def extract_umi_barcode(match_config, barcode_config, reads1, reads2, output_dir, sample, compression_level, correct_barcode, cores=1, batch_size=DEFAULT_BATCH_SIZE, gzip_output=True, cb_len=None, umi_len=None):

    global BARCODE_LEN, UMI_LEN
    if not isinstance(cb_len, int) or isinstance(cb_len, bool) or cb_len <= 0 or cb_len % 2 != 0:
        raise ValueError("cb_len must be a positive even integer")
    if not isinstance(umi_len, int) or isinstance(umi_len, bool) or umi_len <= 0:
        raise ValueError("umi_len must be a positive integer")
    BARCODE_LEN = cb_len // 2
    UMI_LEN = umi_len

    exact_match_stats = [0, 0, 0]
    fuzzy_match_stats = [0, 0, 0]
    n_reads = 0
    n_reads_passed = 0
    overall_start_time = time.time()

    linker1_mm = get_mm_dist(match_config.linker1, match_config.mm_rate)
    linker2_mm = get_mm_dist(match_config.linker2, match_config.mm_rate)
    cores = max(1, int(cores))
    batch_size = max(1, int(batch_size))
    output_suffix = ".fq.gz" if gzip_output else ".fq"

    print("=" * 80)
    print(f"Processing {reads1} and {reads2}")
    print(f"Output files: {output_dir}/{sample}_bc_match_R1{output_suffix} and {output_dir}/{sample}_bc_match_R2{output_suffix}")
    if gzip_output:
        print(f"Compression level: {compression_level}.")
    else:
        print("Output compression: none. This uses more disk space but can speed up preprocessing.")
    print(f"Barcode extraction cores: {cores}.")
    print(f"Barcode extraction batch size: {batch_size}.")
    print(f"Cell barcode length: {cb_len} bp ({BARCODE_LEN} bp per component).")
    print(f"UMI length: {UMI_LEN} bp.")

    if correct_barcode:
        print("Buliding barcode correction maps...")
        bcA_wl = read_whitelist(barcode_config.barcodeA_whitelist)
        bcB_wl = read_whitelist(barcode_config.barcodeB_whitelist)
        print(f"Barcode A whitelist size: {len(bcA_wl)}")
        print(f"Barcode B whitelist size: {len(bcB_wl)}")
        print(f"Total combination {len(bcA_wl) * len(bcB_wl)}")
        barcodeA_correction_map = build_barcode_correction_map(bcA_wl)
        barcodeB_correction_map = build_barcode_correction_map(bcB_wl)
    else:
        print("Skipping barcode correction.")
        barcodeA_correction_map = {}
        barcodeB_correction_map = {}

    worker_context = {
        "match_config": match_config,
        "linker1_mm": linker1_mm,
        "linker2_mm": linker2_mm,
        "correct_barcode": correct_barcode,
        "barcodeA_correction_map": barcodeA_correction_map,
        "barcodeB_correction_map": barcodeB_correction_map,
    }

    open_output = gzip.open if gzip_output else open
    open_kwargs = {"compresslevel": compression_level} if gzip_output else {}

    with tqdm(desc="Extracting barcode/UMI", unit="batch", dynamic_ncols=True) as progress, \
         open_output(f"{output_dir}/{sample}_bc_match_R1{output_suffix}", "wt", **open_kwargs) as out_r1, \
         open_output(f"{output_dir}/{sample}_bc_match_R2{output_suffix}", "wt", **open_kwargs) as out_r2, \
         open_text(reads1) as r1_handle, \
         open_text(reads2) as r2_handle:
        batch_iter = iter_paired_fastq_batches(r1_handle, r2_handle, batch_size)
        if cores == 1:
            for batch in batch_iter:
                result = process_read_batch(batch, worker_context)
                write_batch_result(result, out_r1, out_r2)
                add_batch_stats(result, exact_match_stats, fuzzy_match_stats)
                n_reads += result.n_reads
                n_reads_passed += result.n_reads_passed
                progress.update(1)
                progress.set_postfix(reads=n_reads, passed=n_reads_passed)
        else:
            initargs = (match_config, linker1_mm, linker2_mm, correct_barcode, barcodeA_correction_map, barcodeB_correction_map, BARCODE_LEN, UMI_LEN)
            max_pending = max(cores * 2, 1)
            with ProcessPoolExecutor(max_workers=cores, initializer=init_worker, initargs=initargs) as executor:
                pending = {}
                next_submit_idx = 0
                next_write_idx = 0
                buffered = {}
                exhausted = False

                while pending or not exhausted:
                    while not exhausted and len(pending) < max_pending:
                        try:
                            batch = next(batch_iter)
                        except StopIteration:
                            exhausted = True
                            break
                        future = executor.submit(process_read_batch, batch)
                        pending[future] = next_submit_idx
                        next_submit_idx += 1

                    if not pending:
                        break

                    done, _ = wait(pending, return_when=FIRST_COMPLETED)
                    for future in done:
                        idx = pending.pop(future)
                        buffered[idx] = future.result()

                    while next_write_idx in buffered:
                        result = buffered.pop(next_write_idx)
                        write_batch_result(result, out_r1, out_r2)
                        add_batch_stats(result, exact_match_stats, fuzzy_match_stats)
                        n_reads += result.n_reads
                        n_reads_passed += result.n_reads_passed
                        next_write_idx += 1
                        progress.update(1)
                        progress.set_postfix(reads=n_reads, passed=n_reads_passed)

    print(f"Processed {n_reads} reads, {n_reads_passed} reads passed.")
    percent_passed = n_reads_passed / n_reads * 100 if n_reads else 0
    print(f"Percenatge passed: {percent_passed:.2f}%")
    print(f"Exact match stats: all exact = {exact_match_stats[0]}, linker1 exact = {exact_match_stats[1]}, linker2 exact = {exact_match_stats[2]}")
    print(f"Fuzzy match stats: all fuzzy = {fuzzy_match_stats[0]}, linker1 fuzzy = {fuzzy_match_stats[1]}, linker2 fuzzy = {fuzzy_match_stats[2]}")
    print(f"Overall time: {time.time() - overall_start_time:.2f} seconds.")
    print("=" * 80)


def run_preprocess(
    match_config,
    barcode_config,
    reads1,
    reads2,
    output,
    correct_barcode,
    cb_len,
    umi_len,
    sample,
    compression_level,
    cores,
    batch_size,
    gzip_output,
):
    if cb_len <= 0 or cb_len % 2 != 0:
        raise ValueError("cb_len must be a positive even integer")
    if umi_len <= 0:
        raise ValueError("umi_len must be a positive integer")

    output_path = Path(output) / "fastq_umi_barcode"
    output_path.mkdir(parents=True, exist_ok=True)
    extract_umi_barcode(
        match_config,
        barcode_config,
        reads1,
        reads2,
        output_path,
        sample,
        compression_level,
        correct_barcode,
        cores,
        batch_size,
        gzip_output,
        cb_len,
        umi_len,
    )
    print(output_path)


def main():
    args = parse_args()
    run_preprocess(
        MatchConfig(args.linker1, args.linker2, args.mm_rate),
        BarcodeConfig(
            args.barcodeA_whitelist,
            args.barcodeB_whitelist,
        ),
        args.reads1,
        args.reads2,
        args.output,
        args.correct_barcode,
        args.cb_len,
        args.umi_len,
        args.sample,
        args.compression_level,
        args.core,
        args.batch_size,
        args.gzip_output,
    )


if __name__ == "__main__":
    main()
