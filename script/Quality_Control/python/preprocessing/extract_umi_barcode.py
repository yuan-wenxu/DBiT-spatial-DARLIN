import gzip
from fuzzysearch import find_near_matches
import time
import argparse

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

def find_exact_matches(seq_str, patterns):
    """
    Fast exact matching using Python's built-in string.find()
    Returns: dict with pattern_name -> [positions] mapping
    """
    matches = {}
    for pattern_name, pattern_seq in patterns.items():
        positions, start = [], 0
        while True:
            pos = seq_str.find(pattern_seq, start)
            if pos == -1:
                break
            positions.append(pos)
            start = pos + 1  # allow overlapping
        matches[pattern_name] = positions
    return matches

def find_fuzzy_matches(seq_str, pattern, max_errors):
    """Fallback fuzzy matching (Levenshtein) when exact matching fails."""
    return find_near_matches(pattern, seq_str, max_l_dist=max_errors)

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
    Require exactly one hit for each element.
    """
    result = MatchResult()
    patterns = {'linker1': linker1, 'linker2': linker2}
    exact = find_exact_matches(seq_str, patterns)

    exact_success = (len(exact['linker1']) == 1 and len(exact['linker2']) == 1)
    if exact_success:
        result.linker1_matches = [type('Match', (), {'start': exact['linker1'][0], 'end': exact['linker1'][0] + len(linker1)})()]
        result.linker2_matches = [type('Match', (), {'start': exact['linker2'][0], 'end': exact['linker2'][0] + len(linker2)})()]
        result.match_method = "exact"
        result.match_stats = [1, -1, -1]
        return result

    mixed_success = True
    methods_used = []

    # linker1
    if len(exact['linker1']) == 1:
        result.linker1_matches = [type('Match', (), {'start': exact['linker1'][0], 'end': exact['linker1'][0] + len(linker1)})()]
        methods_used.append('exact')
        result.match_stats[1] = 1
    else:
        result.linker1_matches = find_fuzzy_matches(seq_str, linker1, linker1_mm)
        methods_used.append('fuzzy')
        result.match_stats[1] = 0
        if len(result.linker1_matches) != 1:
            mixed_success = False

    # linker2
    if mixed_success:
        if len(exact['linker2']) == 1:
            result.linker2_matches = [type('Match', (), {'start': exact['linker2'][0], 'end': exact['linker2'][0] + len(linker2)})()]
            methods_used.append('exact')
            result.match_stats[2] = 1
        else:
            result.linker2_matches = find_fuzzy_matches(seq_str, linker2, linker2_mm)
            methods_used.append('fuzzy')
            result.match_stats[2] = 0
            if len(result.linker2_matches) != 1:
                mixed_success = False

    
    if mixed_success:
        result.match_method = "mixed" if 'fuzzy' in methods_used else "exact"
        if result.match_stats[1] == 0 and result.match_stats[2] == 0:
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
    barcodeA = seq[linker2_end: linker2_end + 8]
    barcodeB = seq[linker2_start - 8: linker2_start]
    barcodeA_q = qual[linker2_end: linker2_end + 8]
    barcodeB_q = qual[linker2_start - 8: linker2_start]
    barcode = barcodeB + barcodeA
    barcode_q = barcodeB_q + barcodeA_q
    umi = seq[linker1_end: linker1_end + 10]
    umi_q = qual[linker1_end: linker1_end + 10]
    return barcode, umi, barcode_q, umi_q

def correct_barcode(barcode, barcodeA_correction_map, barcodeB_correction_map):
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

def write_seqrecord_to_fastq(record_id, seq, qual, f):
    """One-shot write of 4 FASTQ lines; quality is expected to be a string."""
    qual_str = qual_to_string(qual)
    f.write(f"@{record_id}\n{seq}\n+\n{qual_str}\n")

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

def main(match_config, barcode_config, reads1, reads2, output_dir, sample, compression_level):

    exact_match_stats = [0, 0, 0, 0]
    fuzzy_match_stats = [0, 0, 0, 0]
    n_reads = 0
    n_reads_passed = 0
    overall_start_time = time.time()

    linker1_mm = get_mm_dist(match_config.linker1, match_config.mm_rate)
    linker2_mm = get_mm_dist(match_config.linker2, match_config.mm_rate)

    print("=" * 80)
    print(f"Processing {reads1} and {reads2}")
    print(f"Output files: {output_dir}/{sample}_bc_match_R1.fq.gz and {output_dir}/{sample}_bc_match_R2.fq.gz")
    print(f"Compression level: {compression_level}.")

    print("Buliding barcode correction maps...")
    bcA_wl = read_whitelist(barcode_config.barcodeA_whitelist)
    bcB_wl = read_whitelist(barcode_config.barcodeB_whitelist)
    print(f"Barcode A whitelist size: {len(bcA_wl)}")
    print(f"Barcode B whitelist size: {len(bcB_wl)}")
    print(f"Total combination {len(bcA_wl) * len(bcB_wl)}")
    barcodeA_correction_map = build_barcode_correction_map(bcA_wl, max_dist = barcode_config.bc_max_dist)
    barcodeB_correction_map = build_barcode_correction_map(bcB_wl, max_dist = barcode_config.bc_max_dist)

    out_r1 = gzip.open(f"{output_dir}/{sample}_bc_match_R1.fq.gz", "wt", compresslevel = compression_level)
    out_r2 = gzip.open(f"{output_dir}/{sample}_bc_match_R2.fq.gz", "wt", compresslevel = compression_level)

    with open_fastq_file(reads1) as r1_handle, open_fastq_file(reads2) as r2_handle:
        for (r1_id, r1_seq, r1_qual), (r2_id, r2_seq, r2_qual) in zip(iter_fastq_raw(r1_handle), iter_fastq_raw(r2_handle)):
            n_reads += 1        
            match_result = find_all_matches(r1_seq, match_config.linker1, match_config.linker2, linker1_mm, linker2_mm)
            if (len(match_result.linker2_matches) == 1) and (len(match_result.linker1_matches) == 1):
                linker2_start = match_result.linker2_matches[0].start
                linker2_end = match_result.linker2_matches[0].end
                linker1_end = match_result.linker1_matches[0].end
                barcode, umi, barcode_q, umi_q = extract_barcode(r1_seq, r1_qual, linker2_start, linker2_end, linker1_end)
                barcode = correct_barcode(barcode, barcodeA_correction_map, barcodeB_correction_map)
                if barcode is None:
                    continue
                if len(barcode) != 16 or len(umi) != 10:
                    continue
                n_reads_passed += 1
                # write matched R1: barcode(16bp)+UMI(10bp)
                new_seq  = barcode + umi
                new_qual = barcode_q + umi_q
                write_seqrecord_to_fastq(r1_id, new_seq, new_qual, out_r1)
                # write matched R2: pass-through
                write_seqrecord_to_fastq(r2_id, r2_seq, r2_qual, out_r2)
            for i in range(3):
                if match_result.match_stats[i] == 0:
                    fuzzy_match_stats[i] += 1
                elif match_result.match_stats[i] == 1:
                    exact_match_stats[i] += 1

    out_r1.close()
    out_r2.close()

    print(f"Processed {n_reads} reads, {n_reads_passed} reads passed.")
    print(f"Percenatge passed: {n_reads_passed / n_reads * 100:.2f}%")
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
    argparser.add_argument('--compression_level', type = int, default = 6, help = 'Compression level for output fastq files')
    argparser.add_argument('--bc_max_dist', type = int, default = 1, help = 'Maximum distance for barcode correction')

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

    match_config = MatchConfig(linker1, linker2, mm_rate)
    barcode = BarcodeConfig(barcodeA_whitelist, barcodeB_whitelist, bc_max_dist)

    main(match_config, barcode, reads1, reads2, output_dir, compression_level)