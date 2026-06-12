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
import math
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


LABELS = ("CA", "RA", "TA")
PLOT_ROOT_NAME = "LR_spatial_plots"
ONLY_OVER_DIRNAME = "only_in_n_LR_gt_count"
BOTH_DIRNAME = "in_both_groups"


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
    parser.add_argument(
        "--cluster-csv",
        required=True,
        help="mRNA cluster CSV with x/y/leiden/color columns used as the background map.",
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
    input_path = input_dir / label / "cellfiltered.csv"
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


def load_cluster_background(cluster_csv: Path) -> pd.DataFrame:
    if not cluster_csv.exists():
        raise SystemExit(f"Cluster CSV not found: {cluster_csv}")

    frame = pd.read_csv(cluster_csv, dtype={"color": str}, keep_default_na=False)
    required_columns = {"x", "y", "leiden", "color"}
    missing_columns = required_columns - set(frame.columns)
    if missing_columns:
        raise SystemExit(
            f"Cluster CSV {cluster_csv} is missing required columns: {missing_columns}"
        )
    frame = sanitize_table(frame)
    frame["color"] = frame["color"].astype(str).str.strip()
    return frame.sort_values(["y", "x"], kind="stable").reset_index(drop=True)


def _hex_to_rgb(color: str) -> tuple[int, int, int]:
    color = color.strip().lstrip("#")
    if len(color) != 6:
        raise ValueError(f"Unsupported color format: {color}")
    return tuple(int(color[i : i + 2], 16) for i in (0, 2, 4))


def choose_contrasting_circle_color(cluster_frame: pd.DataFrame) -> str:
    candidate_colors = [
        "#000000",
        "#ff0055",
        "#00c853",
        "#ffb300",
        "#7c4dff",
        "#00bcd4",
        "#ff6d00",
    ]
    used_colors = []
    for color in cluster_frame["color"].astype(str):
        try:
            used_colors.append(_hex_to_rgb(color))
        except ValueError:
            continue

    if not used_colors:
        return "#000000"

    best_color = "#000000"
    best_score = -1.0
    for candidate in candidate_colors:
        rgb = _hex_to_rgb(candidate)
        min_distance = min(
            math.sqrt(
                (rgb[0] - used[0]) ** 2
                + (rgb[1] - used[1]) ** 2
                + (rgb[2] - used[2]) ** 2
            )
            for used in used_colors
        )
        if min_distance > best_score:
            best_score = min_distance
            best_color = candidate
    return best_color


def _plot_lr_spatial(
    cluster_frame: pd.DataFrame,
    lr_data: pd.DataFrame,
    output_path: Path,
    x_spots: int,
    y_spots: int,
    circle_color: str,
) -> None:
    if {"SR", "LR", "UR"}.issubset(lr_data.columns):
        lr_spots = (
            lr_data.assign(_ur_clean=lr_data["UR"].astype(str).str.strip())
            .groupby(["SR", "LR", "x", "y"], as_index=False)
            .agg(unique_ur_count=("_ur_clean", lambda values: values.replace("", pd.NA).dropna().nunique()))
            .sort_values(["SR", "x", "y"], kind="stable")
            .copy()
        )
    else:
        sort_columns = ["x", "y"]
        if "SR" in lr_data.columns:
            sort_columns = ["SR", "x", "y"]
        lr_spots = (
            lr_data.sort_values(sort_columns, kind="stable")
            .drop_duplicates(subset=["x", "y"], keep="first")
            .copy()
        )
        lr_spots["unique_ur_count"] = 1

    lr_spots["circle_size"] = lr_spots["unique_ur_count"].map(
        lambda count: 30 if count <= 1 else 70 if count == 2 else 140
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(4.8, 4.8))
    ax.scatter(
        cluster_frame["x"],
        cluster_frame["y"],
        s=28,
        c=cluster_frame["color"],
        marker="s",
        alpha=0.95,
        linewidths=0,
    )
    if not lr_spots.empty:
        ax.scatter(
            lr_spots["x"],
            lr_spots["y"],
            s=lr_spots["circle_size"],
            facecolors="none",
            edgecolors=circle_color,
            linewidths=2,
        )

    ax.set_xlim(-0.5, x_spots - 0.5)
    ax.set_ylim(-0.5, y_spots - 0.5)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.tight_layout(pad=0.05)
    fig.savefig(output_path, dpi=160, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def write_lr_plots(
    cluster_frame: pd.DataFrame,
    frame: pd.DataFrame,
    over_limit: pd.DataFrame,
    within_limit: pd.DataFrame,
    output_dir: Path,
    input_stem: str,
) -> tuple[int, int]:
    frame_with_group = []
    if not over_limit.empty:
        part = over_limit.copy()
        part["group"] = "n_LR_gt_count"
        frame_with_group.append(part)
    if not within_limit.empty:
        part = within_limit.copy()
        part["group"] = "n_LR_le_count"
        frame_with_group.append(part)

    if not frame_with_group:
        return 0, 0

    plot_frame = pd.concat(frame_with_group, ignore_index=True)
    over_lrs = set(over_limit["LR"].astype(str).str.strip().replace("", pd.NA).dropna())
    within_lrs = set(within_limit["LR"].astype(str).str.strip().replace("", pd.NA).dropna())
    only_over_lrs = sorted(over_lrs - within_lrs)
    both_lrs = sorted(over_lrs & within_lrs)

    x_spots = max(int(frame["x"].max()), int(cluster_frame["x"].max())) + 1
    y_spots = max(int(frame["y"].max()), int(cluster_frame["y"].max())) + 1
    circle_color = choose_contrasting_circle_color(cluster_frame)

    only_over_dir = output_dir / PLOT_ROOT_NAME / ONLY_OVER_DIRNAME
    both_dir = output_dir / PLOT_ROOT_NAME / BOTH_DIRNAME

    for idx, lr in enumerate(only_over_lrs, start=1):
        lr_data = plot_frame[plot_frame["LR"] == lr].copy()
        sr_count = lr_data["SR"].nunique()
        file_name = f"{input_stem}.only_over.{idx:04d}.sr{sr_count:03d}.png"
        _plot_lr_spatial(
            cluster_frame=cluster_frame,
            lr_data=lr_data,
            output_path=only_over_dir / file_name,
            x_spots=x_spots,
            y_spots=y_spots,
            circle_color=circle_color,
        )

    for idx, lr in enumerate(both_lrs, start=1):
        lr_data = plot_frame[plot_frame["LR"] == lr].copy()
        sr_count = lr_data["SR"].nunique()
        file_name = f"{input_stem}.both_groups.{idx:04d}.sr{sr_count:03d}.png"
        _plot_lr_spatial(
            cluster_frame=cluster_frame,
            lr_data=lr_data,
            output_path=both_dir / file_name,
            x_spots=x_spots,
            y_spots=y_spots,
            circle_color=circle_color,
        )

    return len(only_over_lrs), len(both_lrs)


def process_file(input_path: Path, output_dir: Path, label: str, cluster_frame: pd.DataFrame) -> None:
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

    only_over_lr_count, both_lr_count = write_lr_plots(
        cluster_frame,
        frame,
        over_limit,
        within_limit,
        output_dir,
        input_path.stem,
    )

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
            {
                "label": label,
                "scope": "plot_groups",
                "group": "only_in_n_LR_gt_count",
                "n_rows": int(over_limit[~over_limit["LR"].isin(within_limit["LR"])].shape[0]),
                "n_unique_LR": int(only_over_lr_count),
                "n_unique_SR": int(
                    over_limit[~over_limit["LR"].isin(within_limit["LR"])]["SR"]
                    .astype(str)
                    .str.strip()
                    .replace("", pd.NA)
                    .dropna()
                    .nunique()
                ),
            },
            {
                "label": label,
                "scope": "plot_groups",
                "group": "in_both_groups",
                "n_rows": int(frame[frame["LR"].isin(set(over_limit["LR"]) & set(within_limit["LR"]))].shape[0]),
                "n_unique_LR": int(both_lr_count),
                "n_unique_SR": int(
                    frame[frame["LR"].isin(set(over_limit["LR"]) & set(within_limit["LR"]))]["SR"]
                    .astype(str)
                    .str.strip()
                    .replace("", pd.NA)
                    .dropna()
                    .nunique()
                ),
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
    print(
        f"Saved LR plots: only_in_n_LR_gt_count={only_over_lr_count}, "
        f"in_both_groups={both_lr_count}"
    )


def main() -> None:
    args = parse_args()
    input_dir = Path(args.input_dir)
    output_root = Path(args.output_dir) if args.output_dir else input_dir
    cluster_frame = load_cluster_background(Path(args.cluster_csv))
    output_root.mkdir(parents=True, exist_ok=True)

    for label in args.labels:
        print(f"\n=== Processing {label} ===")
        input_path = resolve_input_path(input_dir, label)
        label_output_dir = output_root / label
        process_file(
            input_path,
            label_output_dir,
            label=label,
            cluster_frame=cluster_frame,
        )

    print("\nDone!")


if __name__ == "__main__":
    main()
