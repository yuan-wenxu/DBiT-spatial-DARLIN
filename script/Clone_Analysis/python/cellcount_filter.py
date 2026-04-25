#!/usr/bin/env python3
"""
Filter DARLIN cellfiltered.csv files so each SR keeps at most `count` unique LR.

For each SR (spot):
- Compute total reads per LR within the spot
- If unique LR <= count, keep all rows
- If unique LR > count, keep only the top `count` LR by total reads

Usage:
  python spot_lr_cellcount_filter.py --input_dir /path/to/data
  python spot_lr_cellcount_filter.py --input_dir /path/to/data --labels RA TA
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


LABELS = ("CA", "RA", "TA")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Filter DARLIN cellfiltered.csv files so the number of unique LR in each "
            "SR does not exceed the predicted cell count."
        )
    )
    parser.add_argument(
        "--input_dir",
        required=True,
        help="Directory containing CA/RA/TA subdirectories.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Directory for outputs (default: input directory).",
    )
    parser.add_argument(
        "--labels",
        nargs="+",
        choices=LABELS,
        default=list(LABELS),
        help="Labels to process.",
    )
    return parser.parse_args()


def filter_spot_lrs_by_cell_count(frame: pd.DataFrame) -> pd.DataFrame:
    """
    For each SR (spot), enforce that the number of unique LR does not exceed
    the predicted cell count (`count`).

    If unique LR > count, keep only the top `count` LR by total reads within
    that SR and drop rows belonging to the rest.
    """
    if frame.empty:
        return frame

    filtered_groups = []

    for _, group in frame.groupby("SR", sort=False):
        cell_count = pd.to_numeric(group["count"], errors="coerce").fillna(0).max()
        cell_count = int(cell_count) if cell_count > 0 else 0

        lr_reads = (
            group.assign(_reads_num=pd.to_numeric(group["reads"], errors="coerce").fillna(0))
            .groupby("LR", as_index=False)["_reads_num"]
            .sum()
            .sort_values("_reads_num", ascending=False)
        )
        unique_lr_count = len(lr_reads)

        if unique_lr_count <= cell_count:
            filtered_groups.append(group)
            continue

        kept_lr = set(lr_reads.head(cell_count)["LR"].tolist())
        filtered_groups.append(group[group["LR"].isin(kept_lr)])

    if not filtered_groups:
        return frame.iloc[0:0].copy()

    return pd.concat(filtered_groups, ignore_index=True)


def process_label(input_dir: Path, output_dir: Path, label: str) -> None:
    input_path = input_dir / label / "cellfiltered.bank_filtered.csv"
    if not input_path.exists():
        raise SystemExit(f"Input file not found: {input_path}")

    frame = pd.read_csv(input_path, dtype=str, keep_default_na=False)
    if frame.empty:
        print(f"{label}: No data")
        return

    required_columns = {"SR", "LR", "reads", "count"}
    missing_columns = required_columns - set(frame.columns)
    if missing_columns:
        raise SystemExit(
            f"Input file {input_path} is missing required columns: {missing_columns}"
        )

    filtered = filter_spot_lrs_by_cell_count(frame)
    print(
        f"{label}: {len(frame)} rows -> {len(filtered)} rows after SR-cellcount LR filtering"
    )

    output_path = output_dir / label / "cellfiltered.bank_filtered.cell_count_filtered.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    filtered.to_csv(output_path, index=False)
    print(f"Saved: {output_path}")


def main() -> None:
    args = parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir) if args.output_dir else input_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    for label in args.labels:
        print(f"\n=== Processing {label} ===")
        process_label(input_dir, output_dir, label)

    print("\nDone!")


if __name__ == "__main__":
    main()