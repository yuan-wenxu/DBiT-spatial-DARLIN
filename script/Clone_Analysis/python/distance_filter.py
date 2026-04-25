#!/usr/bin/env python3
"""
Analyze spatial compactness of clones (LR) from bank-filtered DARLIN data.

For each clone (LR):
1. Collect unique spots (SR) with x/y coordinates
2. Compute all pairwise spot distances
3. Compute mean pairwise distance
4. Classify clone as clustered/dispersed by threshold
5. Plot distance distributions and summary charts

Usage:
  python clone_spatial_distance_analysis.py --input_dir /path/to/L126-S7

Expected structure in input_dir:
  CA/cellfiltered.bank_filtered.csv
  RA/cellfiltered.bank_filtered.csv
  TA/cellfiltered.bank_filtered.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


LABELS = ["CA", "RA", "TA"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute clone-level pairwise spot distances and evaluate whether "
            "clones are spatially clustered or dispersed."
        )
    )
    parser.add_argument(
        "--input_dir",
        required=True,
        help="Path containing CA/RA/TA folders with cellfiltered.bank_filtered.csv.",
    )
    parser.add_argument(
        "--output_dir",
        default=None,
        help="Output directory (default: input_dir).",
    )
    parser.add_argument(
        "--distance-unit",
        choices=["spot", "um"],
        default="um",
        help="Distance unit for reported values and figures.",
    )
    parser.add_argument(
        "--length-spot",
        type=float,
        default=20.0,
        help="Spot size (um).",
    )
    parser.add_argument(
        "--interval",
        type=float,
        default=20.0,
        help="Center-to-center gap contribution (um).",
    )
    parser.add_argument(
        "--cluster-threshold",
        type=float,
        default=200.0,
        help="Threshold for mean pairwise distance to call a clone clustered.",
    )
    parser.add_argument(
    "--labels",
    nargs="+",
    choices=LABELS,
    default=list(LABELS),
    help="Labels to process.",
    )
    return parser.parse_args()


def _pairwise_distances(coords: np.ndarray) -> np.ndarray:
    """Return upper-triangle pairwise Euclidean distances for Nx2 coords."""
    n = coords.shape[0]
    if n < 2:
        return np.array([], dtype=float)

    distances: List[float] = []
    for i in range(n - 1):
        deltas = coords[i + 1 :] - coords[i]
        d = np.sqrt((deltas * deltas).sum(axis=1))
        distances.extend(d.tolist())
    return np.asarray(distances, dtype=float)


def _pairwise_spot_rows(
    lr: str,
    unique_spots: pd.DataFrame,
    scale: float,
) -> list[dict[str, float | int | str]]:
    """Return pairwise distance rows annotated with the contributing SR/x/y pairs."""
    if len(unique_spots) < 2:
        return []

    pair_rows: list[dict[str, float | int | str]] = []
    records = unique_spots[["SR", "x", "y"]].to_dict("records")

    for i in range(len(records) - 1):
        left = records[i]
        left_x = int(left["x"])
        left_y = int(left["y"])

        for j in range(i + 1, len(records)):
            right = records[j]
            right_x = int(right["x"])
            right_y = int(right["y"])
            distance = float(
                np.sqrt((right_x - left_x) ** 2 + (right_y - left_y) ** 2) * scale
            )
            pair_rows.append(
                {
                    "LR": lr,
                    "SR_1": str(left["SR"]),
                    "x_1": left_x,
                    "y_1": left_y,
                    "SR_2": str(right["SR"]),
                    "x_2": right_x,
                    "y_2": right_y,
                    "distance": distance,
                }
            )

    return pair_rows


def _load_filtered_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise SystemExit(f"Input not found: {path}")

    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    need_cols = {"SR", "LR", "x", "y"}
    if not need_cols.issubset(frame.columns):
        missing = need_cols - set(frame.columns)
        raise SystemExit(f"Missing columns in {path}: {missing}")

    frame = frame.copy()
    frame["x"] = pd.to_numeric(frame["x"], errors="coerce")
    frame["y"] = pd.to_numeric(frame["y"], errors="coerce")
    frame = frame.dropna(subset=["x", "y"])
    frame["x"] = frame["x"].astype(int)
    frame["y"] = frame["y"].astype(int)
    return frame


def analyze_label(
    label: str,
    table: pd.DataFrame,
    out_dir: Path,
    args: argparse.Namespace,
) -> Dict[str, int]:
    """Analyze one label and emit CSV/plots."""
    out_dir.mkdir(parents=True, exist_ok=True)
    pitch_um = args.length_spot + args.interval
    scale = pitch_um if args.distance_unit == "um" else 1.0

    clone_rows = []
    pair_rows = []

    grouped = table.groupby("LR", sort=False)
    for lr, group in grouped:
        unique_spots = (
            group[["SR", "x", "y"]]
            .drop_duplicates(subset=["SR", "x", "y"], keep="first")
            .reset_index(drop=True)
        )
        n_spots = len(unique_spots)
        if n_spots < 2:
            clone_rows.append(
                {
                    "LR": lr,
                    "n_spots": n_spots,
                    "pair_count": 0,
                    "mean_pairwise_distance": 0,
                    "median_pairwise_distance": 0,
                    "class": "singleton_or_small",
                }
            )
            continue

        coords = unique_spots[["x", "y"]].to_numpy(dtype=float)
        pair_dist = _pairwise_distances(coords) * scale
        if pair_dist.size == 0:
            clone_rows.append(
                {
                    "LR": lr,
                    "n_spots": n_spots,
                    "pair_count": 0,
                    "mean_pairwise_distance": 0,
                    "median_pairwise_distance": 0,
                    "class": "singleton_or_small",
                }
            )
            continue

        mean_d = float(pair_dist.mean())
        median_d = float(np.median(pair_dist))
        clone_class = "clustered" if mean_d <= args.cluster_threshold else "dispersed"

        clone_rows.append(
            {
                "LR": lr,
                "n_spots": n_spots,
                "pair_count": int(pair_dist.size),
                "mean_pairwise_distance": mean_d,
                "median_pairwise_distance": median_d,
                "class": clone_class,
            }
        )

        pair_rows.extend(_pairwise_spot_rows(lr, unique_spots, scale))

    clone_df = pd.DataFrame(clone_rows)
    pair_df = pd.DataFrame(pair_rows)

    clone_csv = out_dir / f"{label}_clone_distance_summary.csv"
    pair_csv = out_dir / f"{label}_clone_pairwise_distances.csv"
    clone_df.to_csv(clone_csv, index=False)
    pair_df.to_csv(pair_csv, index=False)
    print(f"{label}: saved {clone_csv}")
    print(f"{label}: saved {pair_csv}")

    # Histogram of mean pairwise distance per clone
    valid_clone_dist = clone_df["mean_pairwise_distance"].dropna()
    valid_clone_dist = valid_clone_dist[valid_clone_dist > 0]
    if not valid_clone_dist.empty:
        plt.figure(figsize=(8, 5))
        plt.hist(valid_clone_dist, bins=50, color="#72B7B2", edgecolor="white")
        plt.axvline(
            args.cluster_threshold,
            color="#E45756",
            linestyle="--",
            linewidth=1.5,
            label=f"cluster threshold = {args.cluster_threshold:g}",
        )
        plt.xlabel(f"Mean pairwise distance ({args.distance_unit})")
        plt.ylabel("Number of clones")
        plt.title(f"{label}: Clone mean distance distribution")
        plt.legend()
        plt.tight_layout()
        hist_out = out_dir / f"{label}_clone_mean_distance_distribution.png"
        plt.savefig(hist_out, dpi=300)
        plt.close()
        print(f"{label}: saved {hist_out}")

    # Histogram of all pairwise distances
    if not pair_df.empty:
        plt.figure(figsize=(8, 5))
        plt.hist(pair_df["distance"], bins=50, color="#72B7B2", edgecolor="white")
        plt.axvline(
            args.cluster_threshold,
            color="#E45756",
            linestyle="--",
            linewidth=1.5,
            label=f"cluster threshold = {args.cluster_threshold:g}",
        )
        plt.xlabel(f"Pairwise distance ({args.distance_unit})")
        plt.ylabel("Number of pairs")
        plt.title(f"{label}: All clone pairwise distance distribution")
        plt.legend()
        plt.tight_layout()
        pair_hist_out = out_dir / f"{label}_all_pairwise_distance_distribution.png"
        plt.savefig(pair_hist_out, dpi=300)
        plt.close()
        print(f"{label}: saved {pair_hist_out}")

    # Class summary bar chart
    class_counts = clone_df["class"].value_counts().reindex(
        ["clustered", "dispersed", "singleton_or_small"], fill_value=0
    )
    plt.figure(figsize=(7, 5))
    plt.bar(class_counts.index, class_counts.values, color=["#54A24B", "#EECA3B", "#9D9D9D"])
    plt.ylabel("Number of clones")
    plt.title(f"{label}: Clone compactness classes")
    plt.tight_layout()
    class_out = out_dir / f"{label}_clone_compactness_classes.png"
    plt.savefig(class_out, dpi=300)
    plt.close()
    print(f"{label}: saved {class_out}")

    return {
        "total_clones": int(len(clone_df)),
        "clustered": int(class_counts.get("clustered", 0)),
        "dispersed": int(class_counts.get("dispersed", 0)),
        "small": int(class_counts.get("singleton_or_small", 0)),
    }


def main() -> None:
    args = parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir) if args.output_dir else input_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    all_stats = []
    for label in args.labels:
        input_path = input_dir / label / "cellfiltered.bank_filtered.cell_count_filtered.csv"
        label_out = output_dir / label / "distance_analysis"

        print(f"\n=== {label} ===")
        try:
            table = _load_filtered_table(input_path)
        except SystemExit as exc:
            print(exc)
            continue

        stats = analyze_label(label, table, label_out, args)
        stats["label"] = label
        all_stats.append(stats)

    if all_stats:
        summary = pd.DataFrame(all_stats)[
            ["label", "total_clones", "clustered", "dispersed", "small"]
        ]
        summary_out = output_dir / "clone_distance_summary_all_labels.csv"
        summary.to_csv(summary_out, index=False)
        print(f"\nSaved: {summary_out}")

        print("\nOverall summary:")
        for _, row in summary.iterrows():
            print(
                f"{row['label']}: total={row['total_clones']}, "
                f"clustered={row['clustered']}, dispersed={row['dispersed']}, "
                f"singleton_or_small={row['small']}"
            )
    else:
        print("No labels were successfully analyzed.")


if __name__ == "__main__":
    main()
