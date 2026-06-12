import gzip
from fuzzysearch import find_near_matches
import time
import argparse
from collections import namedtuple
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from itertools import zip_longest

BARCODE_LEN = 8
UMI_LEN = 10
DEFAULT_LINKER_WINDOW = 2
DEFAULT_BATCH_SIZE = 50000
Match = namedtuple("Match", ["start", "end"])
BatchResult = namedtuple("BatchResult", ["r1_records", "r2_records", "n_reads", "n_reads_passed", "exact_match_stats", "fuzzy_match_stats"])

WORKER_CONTEXT = None


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
    # allow 5% mismatch, at least n
    return max(int(len(seq) * rate), n)


def read_whitelist(whitelist_path):
    with open(whitelist_path) as f:
        return set(line.strip() for line in f if line.strip())
    

def build_barcode_correction_map(whitelist, max_dist):
    """
    Build O(1) lookup map for barcode correction.
    For max_dist=1: enumerate all single-substitution neighbors.
    """
    correction_map = {}
    for bc in whitelist:
        correction_map[bc] = bc
    if max_dist >= 1:
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


def open_fastq_file(file_path):
    """Open text or gz FASTQ in text mode."""
    if file_path.endswith(".gz"):
        return gzip.open(file_path, "rt")
    return open(file_path, "r")


def iter_fastq_raw(handle):
    """
    Minimal FASTQ iterator: read 4 lines per record.
    Returns (read_id, seq, qual) as strings. No conversion on quality.
    """
    while True:
        id_line = handle.readline()
        if not id_line:
            break
        seq_line = handle.readline()
        plus_line = handle.readline()
        qual_line = handle.readline()

        if not (seq_line and plus_line and qual_line):
            raise ValueError("Incomplete FASTQ record encountered.")

        if not id_line.startswith("@") or not plus_line.startswith("+"):
            raise ValueError("Invalid FASTQ structure (missing @ or + line).")

        read_id = id_line[1:].strip()
        seq = seq_line.strip()
        qual = qual_line.strip()
        if len(seq) != len(qual):
            raise ValueError(f"Length mismatch (seq {len(seq)} vs qual {len(qual)}) at read {read_id}")
        yield read_id, seq, qual


def iter_paired_fastq_batches(r1_handle, r2_handle, batch_size):
    batch = []
    for r1_record, r2_record in zip_longest(iter_fastq_raw(r1_handle), iter_fastq_raw(r2_handle)):
        if r1_record is None or r2_record is None:
            raise ValueError("R1 and R2 FASTQ files contain different numbers of records.")
        batch.append((r1_record, r2_record))
        if len(batch) >= batch_size:
            yield batch
            batch = []
    if batch:
        yield batch


def get_search_window(seq_len, expected_start, pattern_len, window):
    start = max(0, expected_start - window)
    end = min(seq_len, expected_start + pattern_len + window)
    return start, end


def find_exact_match_near(seq_str, pattern, expected_start, window):
    """Return one exact match in the expected window, or [] if none/ambiguous."""
    window_start, window_end = get_search_window(len(seq_str), expected_start, len(pattern), window)
    sub_seq = seq_str[window_start:window_end]
    first = sub_seq.find(pattern)
    if first == -1:
        return []
    second = sub_seq.find(pattern, first + 1)
    if second != -1:
        return []
    start = window_start + first
    return [Match(start, start + len(pattern))]


def find_fuzzy_matches_near(seq_str, pattern, max_errors, expected_start, window):
    """Fallback fuzzy matching in a small window around the expected linker position."""
    window_start, window_end = get_search_window(len(seq_str), expected_start, len(pattern), window)
    matches = find_near_matches(pattern, seq_str[window_start:window_end], max_l_dist=max_errors)
    return [Match(window_start + match.start, window_start + match.end) for match in matches]


class MatchResult:
    """Container for match results"""
    def __init__(self):
        self.linker1_matches = []
        self.linker2_matches = []
        self.match_method = ""  # "exact", "fuzzy", "mixed", or "failed"
        self.match_stats = [-1, -1, -1]  # [all, linker1, linker2], 0 represents fuzzy, 1 represents exact


def find_all_matches(seq_str, linker1, linker2, linker1_mm, linker2_mm, linker_window=DEFAULT_LINKER_WINDOW):
    """
    Try exact matches first; fall back to fuzzy per element if needed.
    Linker positions are searched only near the expected DBiT amplicon layout:
    barcodeB, linker2, barcodeA, linker1, UMI.
    Require exactly one hit for each element.
    """
    result = MatchResult()

    mixed_success = True
    methods_used = []

    # linker2 starts immediately after barcodeB.
    linker2_expected_start = BARCODE_LEN
    result.linker2_matches = find_exact_match_near(seq_str, linker2, linker2_expected_start, linker_window)
    if len(result.linker2_matches) == 1:
        methods_used.append('exact')
        result.match_stats[2] = 1
    else:
        result.linker2_matches = find_fuzzy_matches_near(seq_str, linker2, linker2_mm, linker2_expected_start, linker_window)
        methods_used.append('fuzzy')
        result.match_stats[2] = 0
        if len(result.linker2_matches) != 1:
            mixed_success = False

    # linker1 starts immediately after barcodeA.
    if mixed_success:
        linker1_expected_start = result.linker2_matches[0].end + BARCODE_LEN
        result.linker1_matches = find_exact_match_near(seq_str, linker1, linker1_expected_start, linker_window)
        if len(result.linker1_matches) == 1:
            methods_used.append('exact')
            result.match_stats[1] = 1
        else:
            result.linker1_matches = find_fuzzy_matches_near(seq_str, linker1, linker1_mm, linker1_expected_start, linker_window)
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
      - barcodeA: 8bp immediately after linker2
      - barcodeB: 8bp immediately before linker2
      - umi: 10bp immediately after linker1
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
    """Use global maps built from whitelists; return corrected 16bp or None."""
    bcB = barcode[:8]
    bcA = barcode[8:]
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
    def __init__(self, barcodeA_whitelist, barcodeB_whitelist, bc_max_dist):
        self.barcodeA_whitelist = barcodeA_whitelist
        self.barcodeB_whitelist = barcodeB_whitelist
        self.bc_max_dist = bc_max_dist


def init_worker(match_config, linker1_mm, linker2_mm, correct_barcode, barcodeA_correction_map, barcodeB_correction_map):
    global WORKER_CONTEXT
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
            if barcode is not None and len(barcode) == 16 and len(umi) == 10:
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


def main(match_config, barcode_config, reads1, reads2, output_dir, sample, compression_level, correct_barcode, cores=1, batch_size=DEFAULT_BATCH_SIZE, gzip_output=True):

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

    if correct_barcode:
        print("Buliding barcode correction maps...")
        bcA_wl = read_whitelist(barcode_config.barcodeA_whitelist)
        bcB_wl = read_whitelist(barcode_config.barcodeB_whitelist)
        print(f"Barcode A whitelist size: {len(bcA_wl)}")
        print(f"Barcode B whitelist size: {len(bcB_wl)}")
        print(f"Total combination {len(bcA_wl) * len(bcB_wl)}")
        barcodeA_correction_map = build_barcode_correction_map(bcA_wl, max_dist = barcode_config.bc_max_dist)
        barcodeB_correction_map = build_barcode_correction_map(bcB_wl, max_dist = barcode_config.bc_max_dist)
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

    with open_output(f"{output_dir}/{sample}_bc_match_R1{output_suffix}", "wt", **open_kwargs) as out_r1, \
         open_output(f"{output_dir}/{sample}_bc_match_R2{output_suffix}", "wt", **open_kwargs) as out_r2, \
         open_fastq_file(reads1) as r1_handle, \
         open_fastq_file(reads2) as r2_handle:
        batch_iter = iter_paired_fastq_batches(r1_handle, r2_handle, batch_size)
        if cores == 1:
            for batch in batch_iter:
                result = process_read_batch(batch, worker_context)
                write_batch_result(result, out_r1, out_r2)
                add_batch_stats(result, exact_match_stats, fuzzy_match_stats)
                n_reads += result.n_reads
                n_reads_passed += result.n_reads_passed
        else:
            initargs = (match_config, linker1_mm, linker2_mm, correct_barcode, barcodeA_correction_map, barcodeB_correction_map)
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

    print(f"Processed {n_reads} reads, {n_reads_passed} reads passed.")
    percent_passed = n_reads_passed / n_reads * 100 if n_reads else 0
    print(f"Percenatge passed: {percent_passed:.2f}%")
    print(f"Exact match stats: all exact = {exact_match_stats[0]}, linker1 exact = {exact_match_stats[1]}, linker2 exact = {exact_match_stats[2]}")
    print(f"Fuzzy match stats: all fuzzy = {fuzzy_match_stats[0]}, linker1 fuzzy = {fuzzy_match_stats[1]}, linker2 fuzzy = {fuzzy_match_stats[2]}")
    print(f"Overall time: {time.time() - overall_start_time:.2f} seconds.")
    print("=" * 80)


if __name__ == '__main__':

    argparser = argparse.ArgumentParser()

    # required arguments
    argparser.add_argument('reads1', type = str, help = 'Path to R1 fastq file')
    argparser.add_argument('reads2', type = str, help = 'Path to R2 fastq file')
    argparser.add_argument('-b1', '--barcodeA_whitelist', type = str, help = 'Path to barcode A whitelist file')
    argparser.add_argument('-b2', '--barcodeB_whitelist', type = str, help = 'Path to barcode B whitelist file')
    argparser.add_argument('-o', '--output', type = str, help = 'Prefix of output fastq files')

    # optional arguments
    argparser.add_argument('-l1', '--linker1', type = str, default = "GTGGCCGATGTTTCGCATCGGCGTACGACT", help = 'Linker 1 sequence')
    argparser.add_argument('-l2', '--linker2', type = str, default = "ATCCACGTGCTTGAGAGGCCAGAGCATTCG", help = 'Linker 2 sequence')
    argparser.add_argument('-m', '--mm_rate', type = float, default = 0.05, help = 'Mismatch rate for linker sequences')
    argparser.add_argument('--compression_level', type = int, default = 6, help = 'Compression level for gzip output files')
    argparser.add_argument('--bc_max_dist', type = int, default = 1, help = 'Maximum distance for barcode correction')
    argparser.add_argument('--correct_barcode', type=str_to_bool, default=False, help='Whether to perform barcode correction')
    argparser.add_argument('--sample', type = str, default = None, help = 'Sample name for output FASTQ files')
    argparser.add_argument('-c', '--cores', type = int, default = 1, help = 'Number of worker processes for barcode extraction')
    argparser.add_argument('--batch_size', type = int, default = DEFAULT_BATCH_SIZE, help = 'Number of read pairs per worker batch')
    argparser.add_argument('--gzip_output', type=str_to_bool, default=True, help='Whether to gzip barcode FASTQ output. False uses more disk space but can be faster.')

    args = argparser.parse_args()

    reads1 = args.reads1  # R1 fastq file
    reads2 = args.reads2  # R2 fastq file
    barcodeA_whitelist = args.barcodeA_whitelist  # barcode A whitelist file
    barcodeB_whitelist = args.barcodeB_whitelist  # barcode B whitelist file
    output_dir = args.output  # output directory

    linker1 = args.linker1  # linker 1 sequence
    linker2 = args.linker2  # linker 2 sequence
    mm_rate = args.mm_rate  # 5% mismatch
    compression_level = args.compression_level  # compression level for output fastq files
    bc_max_dist = args.bc_max_dist  # maximum distance for barcode correction
    correct_barcode = args.correct_barcode  # whether to perform barcode correction
    sample = args.sample
    if sample is None:
        sample = reads1.rsplit("/", 1)[-1].replace("_R1.fq.gz", "").replace("_R1.fastq.gz", "").replace("_R1.fq", "").replace("_R1.fastq", "")

    match_config = MatchConfig(linker1, linker2, mm_rate)
    barcode = BarcodeConfig(barcodeA_whitelist, barcodeB_whitelist, bc_max_dist)

    main(match_config, barcode, reads1, reads2, output_dir, sample, compression_level, correct_barcode, args.cores, args.batch_size, args.gzip_output)
