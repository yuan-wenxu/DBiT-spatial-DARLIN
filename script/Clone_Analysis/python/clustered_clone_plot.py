#!/usr/bin/env python3
"""
Plot all clustered clones using distance_filter outputs.

This script reads:
1. distance_analysis/{label}_clone_distance_summary.csv
2. {label}/cellfiltered.bank_filtered.cell_count_filtered.csv

Then it selects clones with class == "clustered" and renders one spatial plot
for each clustered LR.

Usage:
  python clustered_clone_plot.py --input_dir /path/to/L126-S7
  python clustered_clone_plot.py --input_dir /path/to/L126-S7 --labels RA TA
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image


LABELS = ["CA", "RA", "TA"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot all clustered clones using clone_distance_summary.csv produced "
            "by distance_filter.py."
        )
    )
    parser.add_argument(
        "--input_dir",
        required=True,
        help="Path containing CA/RA/TA folders.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory (default: input_dir).",
    )
    parser.add_argument(
        "--labels",
        nargs="+",
        choices=LABELS,
        default=list(LABELS),
        help="Labels to process.",
    )
    parser.add_argument(
        "--summary-subdir",
        default="distance_analysis",
        help="Subdirectory containing clone_distance_summary.csv files.",
    )
    parser.add_argument(
        "--length-spot",
        type=int,
        default=20,
        help="Practical distance of each spot (um).",
    )
    parser.add_argument(
        "--interval",
        type=int,
        default=20,
        help="Practical interval between spots (um).",
    )
    parser.add_argument(
        "--pixel-length",
        type=float,
        default=0.294,
        help="Pixel size of each spot.",
    )
    parser.add_argument(
        "--x-spots-number",
        type=int,
        default=50,
        help="Number of spots in x direction.",
    )
    parser.add_argument(
        "--y-spots-number",
        type=int,
        default=50,
        help="Number of spots in y direction.",
    )
    return parser.parse_args()


def load_cluster_summary(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise SystemExit(f"Summary file not found: {path}")

    frame = pd.read_csv(path, dtype={"LR": str}, keep_default_na=False)
    required_columns = {"LR", "class", "mean_pairwise_distance", "n_spots"}
    missing_columns = required_columns - set(frame.columns)
    if missing_columns:
        raise SystemExit(f"Missing columns in {path}: {missing_columns}")

    frame = frame.copy()
    frame = frame[frame["class"] == "clustered"].copy()
    frame["mean_pairwise_distance"] = pd.to_numeric(
        frame["mean_pairwise_distance"], errors="coerce"
    )
    frame["n_spots"] = pd.to_numeric(frame["n_spots"], errors="coerce").fillna(0).astype(int)
    frame = frame.sort_values(["mean_pairwise_distance", "n_spots", "LR"], ascending=[True, False, True])
    return frame


def load_filtered_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise SystemExit(f"Input table not found: {path}")

    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    required_columns = {"SR", "LR", "x", "y"}
    missing_columns = required_columns - set(frame.columns)
    if missing_columns:
        raise SystemExit(f"Missing columns in {path}: {missing_columns}")

    frame = frame.copy()
    frame["x"] = pd.to_numeric(frame["x"], errors="coerce")
    frame["y"] = pd.to_numeric(frame["y"], errors="coerce")
    frame = frame.dropna(subset=["x", "y"])
    frame["x"] = frame["x"].astype(int)
    frame["y"] = frame["y"].astype(int)
    return frame


def plot_clone_spatial(
    clone_data: pd.DataFrame,
    label: str,
    clone_name: str,
    output_dir: Path,
    args: argparse.Namespace,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    plot_data = (
        clone_data.sort_values(["SR", "x", "y"])
        .drop_duplicates(subset=["SR", "x", "y"], keep="first")
        .copy()
    )

    x_spots = args.x_spots_number
    y_spots = args.y_spots_number
    length_spot = args.length_spot
    interval = args.interval
    pixel_length = args.pixel_length

    frame_height = int(y_spots * length_spot + (y_spots - 1) * interval)
    frame_width = int(x_spots * length_spot + (x_spots - 1) * interval)
    frame = np.zeros((frame_height, frame_width, 4), dtype=np.uint8)

    for _, row in plot_data.iterrows():
        x_idx = int(row["x"])
        y_idx = int(row["y"])

        x_start = x_idx * (length_spot + interval)
        y_start = y_idx * (length_spot + interval)
        x_end = x_start + length_spot
        y_end = y_start + length_spot

        frame[y_start:y_end, x_start:x_end, 0] = 255
        frame[y_start:y_end, x_start:x_end, 3] = 204

    img = Image.fromarray(frame.astype(np.uint8), mode="RGBA")
    new_width = int(frame_width / pixel_length)
    new_height = int(frame_height / pixel_length)
    img = img.resize((new_width, new_height), resample=Image.NEAREST)

    output_file = output_dir / f"{label}_{clone_name}_spatial.png"
    img.save(output_file)
    print(f"  Plotted: {output_file}")


def process_label(label: str, input_dir: Path, output_dir: Path, args: argparse.Namespace) -> None:
    summary_path = input_dir / label / args.summary_subdir / f"{label}_clone_distance_summary.csv"
    table_path = input_dir / label / "cellfiltered.bank_filtered.cell_count_filtered.csv"

    summary_df = load_cluster_summary(summary_path)
    if summary_df.empty:
        print(f"{label}: No clustered clones found")
        return

    table = load_filtered_table(table_path)
    plot_dir = output_dir / label / args.summary_subdir / "clustered_clone_plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    plotted_rows = []
    for idx, row in enumerate(summary_df.itertuples(index=False), start=1):
        lr = str(row.LR)
        clone_data = table[table["LR"] == lr].copy()
        if clone_data.empty:
            continue

        clone_name = (
            f"cluster_{idx:03d}_spots{int(row.n_spots):02d}_"
            f"mean{float(row.mean_pairwise_distance):.2f}"
        )
        plot_clone_spatial(clone_data, label, clone_name, plot_dir, args)
        plotted_rows.append(
            {
                "LR": lr,
                "n_spots": int(row.n_spots),
                "mean_pairwise_distance": float(row.mean_pairwise_distance),
                "plot_name": clone_name,
            }
        )

    if plotted_rows:
        manifest = pd.DataFrame(plotted_rows)
        manifest_path = plot_dir / f"{label}_clustered_clone_plot_manifest.csv"
        manifest.to_csv(manifest_path, index=False)
        print(f"{label}: saved {manifest_path}")
        print(f"{label}: plotted {len(plotted_rows)} clustered clones")
    else:
        print(f"{label}: No clustered clone plots were generated")


def main() -> None:
    args = parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir) if args.output_dir else input_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    for label in args.labels:
        print(f"\n=== {label} ===")
        try:
            process_label(label, input_dir, output_dir, args)
        except SystemExit as exc:
            print(exc)

    print("\nDone!")


if __name__ == "__main__":
    main()