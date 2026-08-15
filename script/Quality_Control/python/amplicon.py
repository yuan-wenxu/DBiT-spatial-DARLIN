#!/usr/bin/env python3
"""Process extracted DBiT-DARLIN barcodes into a clone table.

Input FASTQs are expected to be paired record-by-record:
  1. SB/UB FASTQ: sequence is SB followed by UB by default.
  2. lineage barcode FASTQ: sequence is the extracted lineage barcode.
"""

import argparse
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm
from umi_tools import UMIClusterer

from utils import (
    add_spatial_coordinates,
    iter_fastq,
    iter_paired_fastq,
    open_text,
    plot_spatial_heatmaps,
    read_barcode_components,
)

PRIMARY_COLOR = "#0072B2"
THRESHOLD_COLOR = "#D55E00"


def str_to_bool(value):
    """Convert common command-line boolean strings to bool."""
    if isinstance(value, bool):
        return value
    if value.lower() in ("yes", "true", "t", "y", "1"):
        return True
    if value.lower() in ("no", "false", "f", "n", "0"):
        return False
    raise argparse.ArgumentTypeError(f"Boolean value expected, got: {value}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert extracted SB/UB and lineage barcode FASTQs to a DARLIN clone table.",
        add_help=False,
    )
    parser.add_argument("--help", action="help", help="Show this help message and exit")
    parser.add_argument(
        "--darlin_reads",
        dest="lineage_bc_fq",
        help="FASTQ(.gz) with DARLIN/lineage barcode sequences.",
    )
    parser.add_argument(
        "--barcode_umi_reads",
        dest="sb_ub_fq",
        required=True,
        help="FASTQ(.gz) with SB+UB sequences.",
    )
    parser.add_argument(
        "--darlin",
        type=str_to_bool,
        default=True,
        help="Whether DARLIN/lineage barcode sequences are available.",
    )
    parser.add_argument(
        "--output_path",
        required=True,
        type=str,
        help="Output directory for final.csv and QC plots.",
    )
    parser.add_argument(
        "--barcodeA_whitelist",
        required=True,
        help="TSV/text file with one barcode A sequence per line.",
    )
    parser.add_argument(
        "--barcodeB_whitelist",
        required=True,
        help="TSV/text file with one barcode B sequence per line.",
    )
    parser.add_argument(
        "--sb-len",
        dest="sb_len",
        type=int,
        required=True,
        help="Length of concatenated spot barcode.",
    )
    parser.add_argument(
        "--ub-len", type=int, required=True, help="Length of UMI barcode."
    )
    parser.add_argument(
        "--x-spots-number",
        "--x_spots_number",
        dest="x_spots_number",
        type=int,
        default=50,
        help="Number of spots in x direction.",
    )
    parser.add_argument(
        "--y-spots-number",
        "--y_spots_number",
        dest="y_spots_number",
        type=int,
        default=50,
        help="Number of spots in y direction.",
    )
    parser.add_argument(
        "--umi_hd_threshold",
        type=int,
        default=1,
        help="Hamming-distance threshold for UMI correction within each SR.",
    )
    parser.add_argument(
        "--min-lb-len",
        type=int,
        default=20,
        help="Minimum lineage barcode length.",
    )
    parser.add_argument(
        "--initial-reads-cutoff",
        type=int,
        default=100,
        help="Minimum reads per raw LB/SB/UB molecule.",
    )
    parser.add_argument(
        "--lb-error-rate",
        dest="lb_error_rate",
        type=float,
        default=0.01,
        help="Lineage barcode correction error rate.",
    )
    parser.add_argument(
        "--lb-min-hd",
        type=int,
        default=1,
        help="Minimum HD threshold for lineage barcode correction.",
    )
    parser.add_argument(
        "--major-fraction-threshold-molecule",
        dest="major_fraction_threshold_molecule",
        type=float,
        default=0.8,
        help="Minimum major LR fraction per SR/UR.",
    )
    parser.add_argument(
        "--reads-fraction-mode",
        dest="reads_fraction_mode",
        choices=("sum", "max"),
        default="sum",
        help="Denominator for major LR filtering within each SR/UR group.",
    )
    parser.add_argument(
        "--slope-cutoff",
        dest="slope_cutoff",
        type=float,
        default=10,
        help="Minimum reads/UMIs per SR.",
    )
    parser.add_argument(
        "--final-reads-cutoff",
        dest="final_reads_cutoff",
        type=int,
        default=10,
        help="Minimum reads per final SR/UR/LR row.",
    )
    args = parser.parse_args()

    if args.darlin and not args.lineage_bc_fq:
        parser.error(
            "amplicon.py requires --darlin_reads when --darlin is True."
        )
    if args.sb_len <= 0 or args.sb_len % 2 != 0:
        parser.error("--sb-len must be a positive even integer.")
    if args.ub_len <= 0:
        parser.error("--ub-len must be a positive integer.")

    args.output_path = Path(args.output_path)
    args.output_path.mkdir(parents=True, exist_ok=True)
    return args


def hamming_dist(first, second):
    if len(first) != len(second):
        return max(len(first), len(second))
    return sum(left != right for left, right in zip(first, second))


def collapse_within_hd(items, max_hd):
    counts = dict(items)
    sequences = sorted(counts, key=lambda sequence: counts[sequence], reverse=True)
    parent = {sequence: sequence for sequence in sequences}
    for index, sequence in enumerate(sequences):
        if parent[sequence] != sequence:
            continue
        for candidate in sequences[index + 1 :]:
            if parent[candidate] != candidate or len(candidate) != len(sequence):
                continue
            if (
                hamming_dist(sequence, candidate) <= max_hd
                and counts[sequence] >= 2 * counts[candidate] - 1
            ):
                parent[candidate] = sequence
    return parent


def build_spot_whitelist(barcode_a_whitelist, barcode_b_whitelist):
    """Build valid concatenated SBs in their observed B+A sequence order."""
    barcode_as = read_barcode_components(barcode_a_whitelist)
    barcode_bs = read_barcode_components(barcode_b_whitelist)
    return {
        barcode_b + barcode_a
        for barcode_b in barcode_bs
        for barcode_a in barcode_as
    }


def neighbors_hd1(sequence):
    return [
        sequence[:index] + base + sequence[index + 1 :]
        for index, current in enumerate(sequence)
        for base in ("A", "C", "G", "T")
        if base != current
    ]


def correct_sb_to_whitelist(data, sb_col="SB", whitelist=None):
    if whitelist is None:
        raise ValueError("whitelist is required for SB correction")
    whitelist = set(whitelist)
    corrected = []
    for barcode in data[sb_col].astype(str):
        if barcode in whitelist:
            corrected.append(barcode)
            continue
        hits = [candidate for candidate in neighbors_hd1(barcode) if candidate in whitelist]
        corrected.append(hits[0] if len(hits) == 1 else None)
    return corrected


def correct_umis(data, sr_col="SR", umi_col="UB", count_col="reads", max_hd=1):
    corrected = []
    clusterer = UMIClusterer(cluster_method="directional")
    for _, subset in tqdm(
        data.groupby(sr_col, sort=False), desc="Correcting UMIs with umi_tools"
    ):
        counts = subset.groupby(umi_col)[count_col].sum().to_dict()
        if not counts:
            continue
        byte_counts = {umi.encode(): int(count) for umi, count in counts.items()}
        mapping = {}
        for group in clusterer(byte_counts, threshold=max_hd):
            representative = max(group, key=lambda umi: byte_counts.get(umi, 0)).decode()
            mapping.update({umi.decode(): representative for umi in group})
        subset = subset.copy()
        subset["UR"] = subset[umi_col].map(lambda umi: mapping.get(umi, umi))
        corrected.append(subset)
    if corrected:
        return pd.concat(corrected, ignore_index=True)
    return data.assign(UR=pd.Series(dtype=str))


def correct_lineage_barcodes(
    data,
    sr_col="SR",
    lb_col="LB",
    lb_len_col="LB_len",
    count_col="reads",
    error_rate=0.01,
    min_hd=1,
):
    data = data.copy()
    if lb_len_col not in data.columns:
        data[lb_len_col] = data[lb_col].astype(str).str.len()
    corrected = []
    groups = data.groupby([sr_col, lb_len_col], sort=False)
    for (_, lb_len), subset in tqdm(groups, desc="Correcting lineage barcodes"):
        counts = subset.groupby(lb_col)[count_col].sum()
        threshold = max(int(round(error_rate * int(lb_len))), min_hd)
        mapping = collapse_within_hd(counts.items(), threshold)
        subset = subset.copy()
        subset["LR"] = subset[lb_col].map(lambda barcode: mapping.get(barcode, barcode))
        corrected.append(subset)
    if corrected:
        return pd.concat(corrected, ignore_index=True)
    return data.assign(LR=pd.Series(dtype=str))


def add_reads_fraction(data, mode):
    data = data.copy()
    if mode == "sum":
        denominator = data.groupby(["SR", "UR"])["reads"].transform("sum")
        data["group_reads"] = denominator
    elif mode == "max":
        denominator = data.groupby(["SR", "UR"])["reads"].transform("max")
        data["group_max_reads"] = denominator
    else:
        raise ValueError(f"Unsupported reads_fraction_mode: {mode}")
    data["reads_fraction"] = data["reads"] / denominator
    return data


def setup_plot_dir(args):
    args.output_path.mkdir(parents=True, exist_ok=True)
    return args.output_path


def finish_plot(figure, output_file):
    for axis in figure.axes:
        axis.xaxis.label.set_size(10)
        axis.yaxis.label.set_size(10)
        axis.tick_params(axis="both", labelsize=10)
        axis.title.set_size(12)
    figure.tight_layout()
    figure.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(figure)
    print(f"Wrote plot: {output_file}")


def plot_lineage_length_hist(lengths, output_file, min_len=None):
    if len(lengths) == 0:
        return
    figure, axis = plt.subplots(figsize=(5, 3))
    axis.hist(
        lengths,
        bins=range(1, 300),
        color=PRIMARY_COLOR,
        edgecolor="white",
        linewidth=0.3,
    )
    if min_len is not None:
        axis.axvline(
            min_len,
            color=THRESHOLD_COLOR,
            linestyle="--",
            linewidth=0.8,
        )
    axis.set(
        xlabel="Sequence Length",
        ylabel="Number of reads",
        title="Distribution of DARLIN Array Sequence\nLengths By Reads",
    )
    finish_plot(figure, output_file)


def get_cutoff_values(data):
    if data.empty:
        return []
    maximum = int(data["reads"].max())
    if maximum < 1:
        return []
    values = list(range(1, min(10, maximum) + 1))
    if maximum >= 11:
        values.extend(range(11, min(50, maximum) + 1, 3))
    if maximum >= 61 and maximum // 2 >= 61:
        values.extend(range(61, maximum // 2 + 1, 10))
    return sorted(set(values))


def plot_reads_cutoff_qc(data, reads_cutoff, output_file):
    cutoffs = get_cutoff_values(data)
    if not cutoffs:
        return
    total_reads = data["reads"].sum()
    molecule_counts = [data.loc[data["reads"] >= cutoff, "UB"].nunique() for cutoff in cutoffs]
    retained = [data.loc[data["reads"] >= cutoff, "reads"].sum() / total_reads for cutoff in cutoffs]
    figure, axes = plt.subplots(2, 1, figsize=(5, 5))
    axes[0].plot(
        cutoffs,
        molecule_counts,
        marker="o",
        markersize=2,
        linewidth=1,
        color=PRIMARY_COLOR,
    )
    axes[0].set(xlabel="Reads Cutoff", ylabel="Number of Molecules", xscale="log", yscale="log")
    axes[1].plot(
        cutoffs,
        retained,
        marker="o",
        markersize=2,
        linewidth=1,
        color=PRIMARY_COLOR,
    )
    axes[1].set(xlabel="Reads Cutoff", ylabel="Frac. of Reads\nRetained", xscale="log", ylim=(0, 1.05))
    for axis in axes:
        axis.axvline(
            reads_cutoff,
            color=THRESHOLD_COLOR,
            linestyle="--",
            linewidth=0.6,
        )
        axis.grid(alpha=0.3)
    finish_plot(figure, output_file)


def plot_reads_fraction_qc(data, threshold, output_file):
    if data.empty:
        return
    figure, axes = plt.subplots(2, 1, figsize=(4, 5))
    axes[0].hist(
        data["reads_fraction"],
        bins=50,
        color=PRIMARY_COLOR,
        edgecolor="white",
        linewidth=0.3,
    )
    axes[0].set(xlabel="Reads Fraction", ylabel="Number of (SR, UR, LR)")
    axes[1].scatter(
        data["reads_fraction"],
        data["reads"],
        s=0.2,
        alpha=0.15,
        color=PRIMARY_COLOR,
    )
    axes[1].set(xlabel="Reads Fraction", ylabel="Reads", yscale="log")
    for axis in axes:
        axis.axvline(
            threshold,
            color=THRESHOLD_COLOR,
            linestyle="--",
            linewidth=0.8,
        )
    finish_plot(figure, output_file)


def plot_sr_reads_umis(summary, output_file):
    if summary.empty:
        return
    figure, axis = plt.subplots(figsize=(4, 3))
    categories = (
        ("<=1", summary["k"] <= 1, "#4575b4"),
        ("<=5", (summary["k"] > 1) & (summary["k"] <= 5), "#91bfdb"),
        ("<=10", (summary["k"] > 5) & (summary["k"] <= 10), "#fee090"),
        (">10", summary["k"] > 10, "#d73027"),
    )
    for _, mask, color in categories:
        subset = summary.loc[mask]
        axis.scatter(subset["n_reads"], subset["n_UR"], s=2, alpha=0.4, color=color)
    axis.set(xscale="log", yscale="log", xlabel="Reads", ylabel="UMIs")
    lower = max(1, summary["n_UR"].min())
    upper = summary["n_UR"].max()
    if lower < upper:
        axis.plot(
            [lower, upper],
            [lower, upper],
            linestyle="--",
            color="red",
            linewidth=1,
        )
    handles = [mpatches.Patch(color=color, label=f"k {label}") for label, _, color in categories]
    axis.legend(handles=handles, title="k = Reads/UMIs", loc="center left", bbox_to_anchor=(1, 0.5), fontsize=8, title_fontsize=9)
    finish_plot(figure, output_file)


def plot_lr_per_sr(data, output_file):
    if data.empty or "n_LR" not in data.columns:
        return
    values = data[["SR", "n_LR"]].drop_duplicates()["n_LR"]
    figure, axis = plt.subplots(figsize=(3, 2))
    axis.hist(
        values,
        bins=range(1, max(8, int(values.max()) + 2)),
        color=PRIMARY_COLOR,
        edgecolor="white",
        linewidth=0.3,
    )
    axis.set(xlabel="Number of LRs per SR", ylabel="Number of SRs", yscale="log")
    finish_plot(figure, output_file)


def read_extracted_fastqs(sb_ub_fq, lineage_bc_fq, sb_len, ub_len):
    """Read amplicon-specific SB/UB and optional lineage FASTQs."""
    if sb_len <= 0 or ub_len <= 0:
        raise ValueError("sb_len and ub_len must be positive integers")
    rows = []
    n_total = 0
    sb_ub_len = sb_len + ub_len

    with open_text(sb_ub_fq) as sb_ub_handle:
        if lineage_bc_fq:
            with open_text(lineage_bc_fq) as lineage_handle:
                records = iter_paired_fastq(sb_ub_handle, lineage_handle)
                for _, sb_ub, _, _, lineage, _ in tqdm(
                    records,
                    desc="Reading extracted FASTQs",
                    unit_scale=True,
                    unit=" reads",
                ):
                    n_total += 1
                    if len(sb_ub) < sb_ub_len:
                        continue
                    rows.append(
                        (
                            lineage,
                            sb_ub[:sb_len],
                            sb_ub[sb_len:sb_ub_len],
                            len(lineage),
                        )
                    )
        else:
            for _, sequence, _ in tqdm(
                iter_fastq(sb_ub_handle),
                desc="Reading extracted SB/UB FASTQ",
                unit_scale=True,
                unit=" reads",
            ):
                n_total += 1
                if len(sequence) < sb_ub_len:
                    continue
                rows.append((sequence[:sb_len], sequence[sb_len:sb_ub_len]))

    print(f"input_reads: {n_total:,}")
    print(f"reads_after_length_filter (barcode+UMI): {len(rows):,}")
    print("\n")
    if lineage_bc_fq:
        return pd.DataFrame(rows, columns=["LB", "SB", "UB", "LB_len"])
    return pd.DataFrame(rows, columns=["SB", "UB"])


def summarize(df, label):
    summary = {
        "reads": int(df["reads"].sum()) if "reads" in df.columns else int(len(df)),
        "molecules": int(len(df)),
    }
    for col in ("SB", "SR", "UB", "UR", "LB", "LR"):
        if col in df.columns:
            print(f"{label}: unique_{col}={int(df[col].nunique()):,}")
    print(f"{label}: reads={summary['reads']:,} molecules={summary['molecules']:,}")
    print("\n")


def apply_sb_correction(df, barcode_a_whitelist, barcode_b_whitelist):
    whitelist = build_spot_whitelist(
        barcode_a_whitelist,
        barcode_b_whitelist,
    )
    observed_whitelist = whitelist & set(df["SB"].unique())
    print(f"Number of whitelist SBs: {len(whitelist):,}")
    print(f"Observed whitelist SBs: {len(observed_whitelist):,}")

    df = df.copy()
    df["SR"] = correct_sb_to_whitelist(df, whitelist=observed_whitelist)
    n_before = len(df)
    df = df[df["SR"].notna()].copy()
    print(f"Rows removed by SB correction: {n_before - len(df):,}")
    summarize(df, "after_sb_correction")
    return df


def process(args):
    plot_dir = setup_plot_dir(args)
    use_lineage = bool(args.darlin)

    df_seq = read_extracted_fastqs(
        args.sb_ub_fq,
        args.lineage_bc_fq if use_lineage else None,
        sb_len=args.sb_len,
        ub_len=args.ub_len,
    )

    if use_lineage:
        plot_lineage_length_hist(df_seq["LB_len"], plot_dir / "lineage_bc_length.png", min_len=args.min_lb_len)
        df = df_seq.groupby(["LB", "SB", "UB", "LB_len"]).size().reset_index(name="reads")
        summarize(df, "collapsed_raw_molecules")
        df = df[df["LB_len"] >= args.min_lb_len].copy()
        summarize_label = "after_length_and_initial_reads_filter"
    else:
        df = df_seq.groupby(["SB", "UB"]).size().reset_index(name="reads")
        summarize(df, "collapsed_raw_molecules")
        summarize_label = "after_initial_reads_filter"

    plot_reads_cutoff_qc(df, args.initial_reads_cutoff, plot_dir / "reads_cutoff_qc.png")
    df = df[df["reads"] >= args.initial_reads_cutoff].copy()
    df.sort_values(by="reads", ascending=False, inplace=True)
    summarize(df, summarize_label)

    df = apply_sb_correction(
        df,
        args.barcodeA_whitelist,
        args.barcodeB_whitelist,
    )
    df = correct_umis(df, max_hd=args.umi_hd_threshold)

    if use_lineage:
        df = (
            df.groupby(["LB", "SR", "UR"], as_index=False)
            .agg(reads=("reads", "sum"))
        )
        df["LB_len"] = df["LB"].str.len()
    else:
        df = df.groupby(["SR", "UR"], as_index=False).agg(reads=("reads", "sum"))

    df.sort_values(by="reads", ascending=False, inplace=True)
    summarize(df, "after_umi_correction")

    if use_lineage:
        df = correct_lineage_barcodes(
            df,
            error_rate=args.lb_error_rate,
            min_hd=args.lb_min_hd,
        )
        df = df.groupby(["SR", "UR", "LR"], as_index=False).agg(reads=("reads", "sum"))
        df.sort_values(by="reads", ascending=False, inplace=True)
        summarize(df, "after_lineage_correction")

        df = add_reads_fraction(df, args.reads_fraction_mode)
        plot_reads_fraction_qc(df, args.major_fraction_threshold_molecule, plot_dir / "reads_fraction_qc.png")
        df_major = df[df["reads_fraction"] >= args.major_fraction_threshold_molecule].copy()
        reads_removed_as_amplification_error = int(df["reads"].sum() - df_major["reads"].sum())
        print(f"reads_removed_as_amplification_error: {reads_removed_as_amplification_error:,}")
        summarize(df_major, "after_major_lr_filter")
    else:
        df_major = df

    sr_summary = df_major.groupby("SR").agg(n_reads=("reads", "sum"), n_UR=("UR", "nunique")).reset_index()
    sr_summary["k"] = sr_summary["n_reads"] / sr_summary["n_UR"]
    plot_sr_reads_umis(sr_summary, plot_dir / "sr_reads_vs_umis.png")

    df_final = df_major.merge(sr_summary[["SR", "k"]], on="SR", how="left")
    df_final = df_final[(df_final["k"] >= args.slope_cutoff) & (df_final["reads"] >= args.final_reads_cutoff)].copy()
    reads_removed_as_capture_oligo_carryover = int(df_major["reads"].sum() - df_final["reads"].sum())
    print(f"reads_removed_as_capture_oligo_carryover: {reads_removed_as_capture_oligo_carryover:,}")
    summarize(df_final, "after_low_quality_sr_filter")
    if use_lineage:
        df_final["n_LR"] = df_final.groupby("SR")["LR"].transform("nunique")
        plot_lr_per_sr(df_final, plot_dir / "lr_per_sr_hist.png")

    df_final = add_spatial_coordinates(
        df_final,
        "SR",
        args.barcodeA_whitelist,
        args.sb_len,
        barcode_b_whitelist_path=args.barcodeB_whitelist,
    )
    out_final = args.output_path / "final.csv"
    df_final.to_csv(out_final, index=False)
    print(f"Wrote final table: {out_final.resolve()}")
    plot_spatial_heatmaps(
        out_final,
        args.barcodeA_whitelist,
        args.barcodeB_whitelist,
        args.output_path,
        args.sb_len,
        args.x_spots_number,
        args.y_spots_number,
    )


def main():
    args = parse_args()
    process(args)


if __name__ == "__main__":
    main()
