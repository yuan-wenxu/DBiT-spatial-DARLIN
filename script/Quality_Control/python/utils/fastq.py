import gzip

import pandas as pd
from tqdm import tqdm


def open_fastq_file(file_path):
    """Open text or gz FASTQ/text files in text mode."""
    file_path = str(file_path)
    if file_path.endswith(".gz"):
        return gzip.open(file_path, "rt")
    return open(file_path, "r")


def iter_fastq_paired(handle1, handle2):
    """Yield paired FASTQ records from two already-open handles."""
    while True:
        id_line1 = handle1.readline()
        if not id_line1:
            break
        seq_line1 = handle1.readline()
        plus_line1 = handle1.readline()
        qual_line1 = handle1.readline()

        if not (seq_line1 and plus_line1 and qual_line1):
            raise ValueError("Incomplete FASTQ record encountered in file 1.")
        if not id_line1.startswith("@") or not plus_line1.startswith("+"):
            raise ValueError("Invalid FASTQ structure (missing @ or + line) in file 1.")

        id_line2 = handle2.readline()
        if not id_line2:
            raise ValueError("File 2 ended before file 1.")
        seq_line2 = handle2.readline()
        plus_line2 = handle2.readline()
        qual_line2 = handle2.readline()

        if not (seq_line2 and plus_line2 and qual_line2):
            raise ValueError("Incomplete FASTQ record encountered in file 2.")
        if not id_line2.startswith("@") or not plus_line2.startswith("+"):
            raise ValueError("Invalid FASTQ structure (missing @ or + line) in file 2.")

        read_id1 = id_line1[1:].strip()
        seq1 = seq_line1.strip()
        qual1 = qual_line1.strip()
        if len(seq1) != len(qual1):
            raise ValueError(f"Length mismatch (seq {len(seq1)} vs qual {len(qual1)}) at read {read_id1} in file 1")

        read_id2 = id_line2[1:].strip()
        seq2 = seq_line2.strip()
        qual2 = qual_line2.strip()
        if len(seq2) != len(qual2):
            raise ValueError(f"Length mismatch (seq {len(seq2)} vs qual {len(qual2)}) at read {read_id2} in file 2")

        yield read_id1, seq1, qual1, read_id2, seq2, qual2


def read_extracted_fastqs(sb_ub_fq, lineage_bc_fq, sb_len=16, ub_len=10):
    rows = []
    n_total = 0
    n_kept = 0
    sb_ub_len = sb_len + ub_len

    with open_fastq_file(sb_ub_fq) as fq1, open_fastq_file(lineage_bc_fq) as fq2:
        iterator = iter_fastq_paired(fq1, fq2)
        for _read_id1, seq1, _qual1, _read_id2, seq2, _qual2 in tqdm(
            iterator,
            desc="Reading extracted FASTQs",
            unit_scale=True,
            unit=" reads",
        ):
            n_total += 1
            if len(seq1) < sb_ub_len:
                continue
            sb = seq1[:sb_len]
            ub = seq1[sb_len:sb_ub_len]
            lb = seq2
            rows.append((lb, sb, ub, len(lb)))
            n_kept += 1

    print(f"input_reads: {n_total:,}")
    print(f"reads_after_length_filter (barcode+UMI): {n_kept:,}")
    print("\n")
    return pd.DataFrame(rows, columns=["LB", "SB", "UB", "LB_len"])
