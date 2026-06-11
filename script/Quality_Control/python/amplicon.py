#!/usr/bin/env python3
"""Process extracted DBiT-DARLIN barcodes into a clone table.

Input FASTQs are expected to be paired record-by-record:
  1. SB/UB FASTQ: sequence is SB followed by UB by default.
  2. lineage barcode FASTQ: sequence is the extracted lineage barcode.
"""

import argparse
from pathlib import Path

from plot import (
    plot_sr_reads_umis,
    plot_lineage_length_hist,
    plot_lr_per_sr,
    plot_reads_cutoff_qc,
    plot_reads_fraction_qc,
    setup_plot_dir,
)
from utils import (
    add_reads_fraction,
    correct_sb_to_whitelist,
    correct_lineage_barcodes,
    correct_umis,
    read_extracted_fastqs,
    read_whitelist,
)


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
        description="Convert extracted SB/UB and lineage barcode FASTQs to a DARLIN clone table."
    )
    parser.add_argument(
        "-dr",
        "--darlin_reads",
        dest="lineage_bc_fq",
        required=True,
        help="FASTQ(.gz) with DARLIN/lineage barcode sequences.",
    )
    parser.add_argument(
        "-bu",
        "--barcode_umi_reads",
        "--bc_umi_reads",
        dest="sb_ub_fq",
        required=True,
        help="FASTQ(.gz) with SB+UB sequences.",
    )
    parser.add_argument(
        "-d",
        "--darlin",
        type=str_to_bool,
        default=True,
        help="Whether DARLIN sequences are provided. This script requires True.",
    )
    parser.add_argument(
        "-o",
        "--output_path",
        required=True,
        type=str,
        help="Output directory for final.csv and QC plots.",
    )
    parser.add_argument("--whitelist", help="Optional TSV/text file with one DBiT barcode per line for SB correction.")
    parser.add_argument("--sb-len", dest="sb_len", type=int, default=16, help="Length of concatenated spot barcode.")
    parser.add_argument("--ub-len", type=int, default=10, help="Length of UMI barcode.")
    parser.add_argument("--umi_hd_threshold", type=int, default=1, help="Hamming-distance threshold for UMI correction within each SR.")
    parser.add_argument("--min-lb-len", type=int, default=20, help="Minimum lineage barcode length.")
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
    parser.add_argument("--lb-min-hd", type=int, default=1, help="Minimum HD threshold for lineage barcode correction.")
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

    if not args.darlin:
        parser.error("amplicon.py requires -d/--darlin True because lineage barcode FASTQ is required.")

    args.output_path = Path(args.output_path)
    args.output_path.mkdir(parents=True, exist_ok=True)

    return args


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


def apply_sb_correction(df, whitelist):
    if whitelist is None:
        raise ValueError("whitelist is required for SB correction")
    if isinstance(whitelist, (str, Path)):
        whitelist = read_whitelist(whitelist)
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

    df_seq = read_extracted_fastqs(
        args.sb_ub_fq,
        args.lineage_bc_fq,
        sb_len=args.sb_len,
        ub_len=args.ub_len,
    )

    plot_lineage_length_hist(df_seq["LB_len"], plot_dir / "lineage_bc_length.png", min_len=args.min_lb_len)

    df = df_seq.groupby(["LB", "SB", "UB", "LB_len"]).size().reset_index(name="reads")
    summarize(df, "collapsed_raw_molecules")

    df = df[df["LB_len"] >= args.min_lb_len].copy()
    plot_reads_cutoff_qc(df, args.initial_reads_cutoff, plot_dir / "reads_cutoff_qc.png")
    df = df[df["reads"] >= args.initial_reads_cutoff].copy()
    df.sort_values(by="reads", ascending=False, inplace=True)
    summarize(df, "after_length_and_initial_reads_filter")

    df = apply_sb_correction(df, args.whitelist)
    df = correct_umis(df, max_hd=args.umi_hd_threshold)
    df = (
        df.groupby(["LB", "SR", "UR"], as_index=False)
        .agg(reads=("reads", "sum"))
    )
    df["LB_len"] = df["LB"].str.len()
    df.sort_values(by="reads", ascending=False, inplace=True)
    summarize(df, "after_umi_correction")

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

    sr_summary = df_major.groupby("SR").agg(n_reads=("reads", "sum"), n_UR=("UR", "nunique")).reset_index()
    sr_summary["k"] = sr_summary["n_reads"] / sr_summary["n_UR"]
    plot_sr_reads_umis(sr_summary, plot_dir / "sr_reads_vs_umis.png")

    df_final = df_major.merge(sr_summary[["SR", "k"]], on="SR", how="left")
    df_final = df_final[(df_final["k"] >= args.slope_cutoff) & (df_final["reads"] >= args.final_reads_cutoff)].copy()
    df_final["n_LR"] = df_final.groupby("SR")["LR"].transform("nunique")
    reads_removed_as_capture_oligo_carryover = int(df_major["reads"].sum() - df_final["reads"].sum())
    print(f"reads_removed_as_capture_oligo_carryover: {reads_removed_as_capture_oligo_carryover:,}")
    summarize(df_final, "after_low_quality_sr_filter")
    plot_lr_per_sr(df_final, plot_dir / "lr_per_sr_hist.png")

    out_final = args.output_path / "final.csv"
    df_final.to_csv(out_final, index=False)
    print(f"Wrote final table: {out_final.resolve()}")


def main():
    args = parse_args()
    process(args)


if __name__ == "__main__":
    main()
