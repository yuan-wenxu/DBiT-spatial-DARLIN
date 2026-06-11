import pandas as pd
from tqdm import tqdm
from umi_tools import UMIClusterer

from .fastq import open_fastq_file


def hamming_dist(a, b):
    """Return Hamming distance; unequal lengths are treated as very distant."""
    if len(a) != len(b):
        return max(len(a), len(b))
    return sum(c1 != c2 for c1, c2 in zip(a, b))


def collapse_within_hd(items, max_hd):
    """Directional collapse of sequences within a Hamming-distance threshold."""
    counts = dict(items)
    seqs = sorted(counts, key=lambda seq: counts[seq], reverse=True)
    parent = {seq: seq for seq in seqs}

    for i, seq in enumerate(seqs):
        if parent[seq] != seq:
            continue
        count_hi = counts[seq]
        for candidate in seqs[i + 1:]:
            if parent[candidate] != candidate:
                continue
            if len(candidate) != len(seq):
                continue
            if hamming_dist(seq, candidate) <= max_hd:
                count_lo = counts[candidate]
                if count_hi >= 2 * count_lo - 1:
                    parent[candidate] = seq
    return parent


def read_whitelist(whitelist_file):
    with open_fastq_file(whitelist_file) as handle:
        barcodes = [line.strip() for line in handle if line.strip()]
    return set(i + j for i in barcodes for j in barcodes)


def neighbors_hd1(seq):
    bases = ("A", "C", "G", "T")
    out = []
    for i, ch in enumerate(seq):
        for base in bases:
            if base != ch:
                out.append(seq[:i] + base + seq[i + 1:])
    return out


def correct_sb_to_whitelist(df, sb_col="SB", whitelist=None):
    if whitelist is None:
        raise ValueError("whitelist is required for SB correction")
    whitelist = set(whitelist)
    corrected = []

    for sb in df[sb_col].astype(str):
        if sb in whitelist:
            corrected.append(sb)
            continue
        hits = [candidate for candidate in neighbors_hd1(sb) if candidate in whitelist]
        corrected.append(hits[0] if len(hits) == 1 else None)

    return corrected


def correct_umis(df, sr_col="SR", umi_col="UB", count_col="reads", max_hd=1):
    df = df.copy()
    corrected = []
    clusterer = UMIClusterer(cluster_method="directional")

    for _sr, sub in tqdm(df.groupby(sr_col, sort=False), desc="Correcting UMIs with umi_tools"):
        umi_counts = sub.groupby(umi_col)[count_col].sum().to_dict()
        if not umi_counts:
            continue

        umi_counts_bytes = {umi.encode(): int(count) for umi, count in umi_counts.items()}
        umi_groups = clusterer(umi_counts_bytes, threshold=max_hd)
        umi_to_rep = {}
        for group in umi_groups:
            representative = max(group, key=lambda umi: umi_counts_bytes.get(umi, 0)).decode()
            for umi in group:
                umi_to_rep[umi.decode()] = representative

        sub = sub.copy()
        sub["UR"] = sub[umi_col].map(lambda umi: umi_to_rep.get(umi, umi))
        corrected.append(sub)

    return pd.concat(corrected, ignore_index=True) if corrected else df.assign(UR=pd.Series(dtype=str))


def correct_lineage_barcodes(
    df,
    sr_col="SR",
    lb_col="LB",
    lb_len_col="LB_len",
    count_col="reads",
    error_rate=0.01,
    min_hd=1,
):
    df = df.copy()
    if lb_len_col not in df.columns:
        df[lb_len_col] = df[lb_col].astype(str).str.len()

    corrected = []
    for (_sr, lb_len), sub in tqdm(
        df.groupby([sr_col, lb_len_col], sort=False),
        desc="Correcting lineage barcodes",
    ):
        counts = sub.groupby(lb_col)[count_col].sum()
        hd_threshold = max(int(round(error_rate * int(lb_len))), min_hd)
        lb_to_rep = collapse_within_hd(counts.items(), max_hd=hd_threshold)

        sub = sub.copy()
        sub["LR"] = sub[lb_col].map(lambda lb: lb_to_rep.get(lb, lb))
        corrected.append(sub)

    return pd.concat(corrected, ignore_index=True) if corrected else df.assign(LR=pd.Series(dtype=str))


def add_reads_fraction(df, mode):
    df = df.copy()
    group_cols = ["SR", "UR"]

    if mode == "sum":
        df["group_reads"] = df.groupby(group_cols)["reads"].transform("sum")
        df["reads_fraction"] = df["reads"] / df["group_reads"]
    elif mode == "max":
        df["group_max_reads"] = df.groupby(group_cols)["reads"].transform("max")
        df["reads_fraction"] = df["reads"] / df["group_max_reads"]
    else:
        raise ValueError(f"Unsupported reads_fraction_mode: {mode}")

    return df
