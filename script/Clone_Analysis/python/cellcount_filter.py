#!/usr/bin/env python3
"""
Split DARLIN filtered CSV rows by whether each SR exceeds its predicted cell count.

For each SR (spot):
- Compute the number of unique LR within the spot
- Compare that number with the predicted cell count (`count`)
- Write one table for SR where n_LR > count
- Write another table for SR where n_LR <= count

Usage:
  python cellcount_filter.py --input-dir /path/to/L126-S9 --labels RA
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


LABELS = ("CA", "RA", "TA")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Split DARLIN filtered CSV rows according to whether each SR has more "
            "unique LR than the predicted cell count."
        )
    )
    parser.add_argument(
        "--input-dir",
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


def classify_spots_by_cell_count(
    frame: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Split rows into two tables according to whether the number of unique LR in each
    SR exceeds the predicted cell count, and return SR-level summary statistics.
    """
    if frame.empty:
        empty = frame.iloc[0:0].copy()
        return empty, empty, pd.DataFrame(
            columns=["SR", "cell_count", "unique_LR_count", "row_count", "extra_LR_count", "group"]
        )

    summary_rows = []
    for sr, group in frame.groupby("SR", sort=False):
        cell_count = pd.to_numeric(group["count"], errors="coerce").fillna(0).max()
        cell_count = int(cell_count) if cell_count > 0 else 0

        unique_lr_count = (
            group["LR"]
            .astype(str)
            .str.strip()
            .replace("", pd.NA)
            .dropna()
            .nunique()
        )
        summary_rows.append(
            {
                "SR": sr,
                "cell_count": cell_count,
                "unique_LR_count": int(unique_lr_count),
                "row_count": int(len(group)),
                "extra_LR_count": int(max(unique_lr_count - cell_count, 0)),
                "group": "n_LR_gt_count" if unique_lr_count > cell_count else "n_LR_le_count",
            }
        )

    summary = pd.DataFrame(summary_rows)
    annotated = frame.merge(summary, on="SR", how="left")

    over_limit = (
        annotated[annotated["group"] == "n_LR_gt_count"]
        .drop(columns=["cell_count", "unique_LR_count", "row_count", "extra_LR_count", "group"])
        .copy()
    )
    within_limit = (
        annotated[annotated["group"] == "n_LR_le_count"]
        .drop(columns=["cell_count", "unique_LR_count", "row_count", "extra_LR_count", "group"])
        .copy()
    )
    summary = summary.sort_values(["SR"], kind="stable").reset_index(drop=True)
    return over_limit, within_limit, summary


def resolve_input_path(input_dir: Path, label: str) -> Path:
    input_path = input_dir / label / "tissuefiltered.csv"
    if input_path.exists():
        return input_path
    raise SystemExit(f"Input file not found: {input_path}")


def sanitize_table(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    frame["x"] = pd.to_numeric(frame["x"], errors="coerce")
    frame["y"] = pd.to_numeric(frame["y"], errors="coerce")
    frame = frame.dropna(subset=["x", "y"])
    frame["x"] = frame["x"].astype(int)
    frame["y"] = frame["y"].astype(int)
    return frame


def process_file(input_path: Path, output_dir: Path, label: str) -> None:
    if not input_path.exists():
        raise SystemExit(f"Input file not found: {input_path}")

    frame = pd.read_csv(input_path, dtype=str, keep_default_na=False)
    if frame.empty:
        print(f"{label}: No data")
        return

    required_columns = {"SR", "UR", "LR", "reads", "count", "x", "y"}
    missing_columns = required_columns - set(frame.columns)
    if missing_columns:
        raise SystemExit(
            f"Input file {input_path} is missing required columns: {missing_columns}"
        )

    frame = sanitize_table(frame)
    over_limit, within_limit, sr_summary = classify_spots_by_cell_count(frame)
    over_limit = over_limit.sort_values(["SR"], kind="stable").reset_index(drop=True)
    within_limit = within_limit.sort_values(["SR"], kind="stable").reset_index(drop=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    over_path = output_dir / f"{input_path.stem}.n_LR_gt_count.csv"
    within_path = output_dir / f"{input_path.stem}.n_LR_le_count.csv"
    summary_path = output_dir / f"{input_path.stem}.count_summary.txt"
    over_limit.to_csv(over_path, index=False)
    within_limit.to_csv(within_path, index=False)

    overall_summary = pd.DataFrame(
        [
            {
                "label": label,
                "scope": "overall",
                "group": "all_rows",
                "n_rows": int(len(frame)),
                "n_unique_LR": int(frame["LR"].astype(str).str.strip().replace("", pd.NA).dropna().nunique()),
                "n_unique_SR": int(frame["SR"].astype(str).str.strip().replace("", pd.NA).dropna().nunique()),
            },
            {
                "label": label,
                "scope": "overall",
                "group": "n_LR_le_count",
                "n_rows": int(len(within_limit)),
                "n_unique_LR": int(within_limit["LR"].astype(str).str.strip().replace("", pd.NA).dropna().nunique()),
                "n_unique_SR": int(within_limit["SR"].astype(str).str.strip().replace("", pd.NA).dropna().nunique()),
            },
            {
                "label": label,
                "scope": "overall",
                "group": "n_LR_gt_count",
                "n_rows": int(len(over_limit)),
                "n_unique_LR": int(over_limit["LR"].astype(str).str.strip().replace("", pd.NA).dropna().nunique()),
                "n_unique_SR": int(over_limit["SR"].astype(str).str.strip().replace("", pd.NA).dropna().nunique()),
            },
        ]
    )

    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write("# Overall summary\n")
        handle.write(overall_summary.to_csv(index=False))
        handle.write("\n")
        handle.write("# SR details where n_LR > count\n")
        if sr_summary[sr_summary["group"] == "n_LR_gt_count"].empty:
            handle.write("None\n")
        else:
            handle.write(
                sr_summary[sr_summary["group"] == "n_LR_gt_count"].to_csv(index=False)
            )

    print(f"Saved: {over_path}")
    print(f"Saved: {within_path}")
    print(f"Saved: {summary_path}")


def main() -> None:
    args = parse_args()
    input_dir = Path(args.input_dir)
    output_root = Path(args.output_dir) if args.output_dir else input_dir
    output_root.mkdir(parents=True, exist_ok=True)

    for label in args.labels:
        print(f"\n=== Processing {label} ===")
        input_path = resolve_input_path(input_dir, label)
        label_output_dir = output_root / label
        process_file(
            input_path,
            label_output_dir,
            label=label,
        )

    print("\nDone!")


if __name__ == "__main__":
    main()
