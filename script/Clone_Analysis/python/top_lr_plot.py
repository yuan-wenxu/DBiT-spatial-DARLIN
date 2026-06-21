#!/usr/bin/env python3
"""
Plot the largest LR entries from allele-bank-filtered CSV files.

"Largest" is defined as the LR observed in the greatest number of unique SR
(spots). The plotting style matches the current Clone_Analysis overlays:
- cluster background from an mRNA data_cellfiltered.csv
- a contrasting hollow circle on spots containing the selected LR
- circle size split into 3 bins by unique UR count per SR/LR

Usage:
  python top_lr_plot.py --input-dir /path/to/L126-S7 --labels RA --top-n 10
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import textwrap

import anndata as ad
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import pandas as pd


LABELS = ("CA", "RA", "TA")
PLOT_EDGE_PAD = 0.9
TITLE_WRAP_WIDTH = 36
TITLE_LINE_COUNT = 3
FIG_SIZE = (5.0, 5.1)
FIG_TOP = 0.84
FIG_BOTTOM = 0.04


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot the largest LR entries from allele-bank-filtered CSV files."
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
        "--top-n",
        type=int,
        default=10,
        help="Number of LR plots to generate per label.",
    )
    parser.add_argument(
        "--cluster-csv",
        required=True,
        help=(
            "mRNA cluster CSV with x/y/leiden/color columns."
        ),
    )
    parser.add_argument(
        "--rotate",
        type=int,
        choices=[0, 90, 180, 270],
        default=0,
        help="Rotate the plot clockwise by 0, 90, 180, or 270 degrees.",
    )
    parser.add_argument(
        "--cluster-alpha",
        type=float,
        default=0.5,
        help="Alpha for the Leiden cluster background spots. Default: 0.5.",
    )
    return parser.parse_args()


def resolve_input_path(input_dir: Path, label: str) -> Path:
    input_path = input_dir / label / "cellfiltered.bank_filtered.csv"
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
    required_columns = {"x", "y", "leiden"}
    missing_columns = required_columns - set(frame.columns)
    if missing_columns:
        raise SystemExit(
            f"Cluster CSV {cluster_csv} is missing required columns: {missing_columns}"
        )
    if "color" not in frame.columns:
        h5ad_path = cluster_csv.parent / "clustered.h5ad"
        leiden_values = sorted(
            pd.Series(frame["leiden"]).astype(str).drop_duplicates().tolist(),
            key=lambda value: int(value) if value.isdigit() else value,
        )
        color_map: dict[str, str] = {}
        if h5ad_path.exists():
            adata = ad.read_h5ad(h5ad_path, backed="r")
            if "leiden_colors" in adata.uns:
                leiden_colors = list(adata.uns["leiden_colors"])
                color_map = {
                    leiden: str(leiden_colors[idx])
                    for idx, leiden in enumerate(leiden_values)
                    if idx < len(leiden_colors)
                }

        if not color_map:
            fallback_colors = [mcolors.to_hex(color) for color in plt.cm.tab20.colors]
            color_map = {
                leiden: fallback_colors[idx % len(fallback_colors)]
                for idx, leiden in enumerate(leiden_values)
            }

        frame["color"] = frame["leiden"].astype(str).map(color_map).fillna("#bdbdbd")

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


def legend_marker_size_from_scatter_area(area: float) -> float:
    return math.sqrt(area) * 0.95


def format_fixed_height_title(text: str, width: int = TITLE_WRAP_WIDTH, line_count: int = TITLE_LINE_COUNT) -> str:
    lines = textwrap.wrap(text, width=width) or [""]
    if len(lines) > line_count:
        lines = lines[:line_count]
        lines[-1] = lines[-1][: max(width - 3, 0)].rstrip() + "..."
    lines.extend([" "] * (line_count - len(lines)))
    return "\n".join(lines)


def rotate_coordinates(
    frame: pd.DataFrame,
    x_spots: int,
    y_spots: int,
    angle: int,
) -> tuple[pd.DataFrame, int, int]:
    rotated = frame.copy()
    if angle == 0:
        rotated["plot_x"] = rotated["x"]
        rotated["plot_y"] = rotated["y"]
        return rotated, x_spots, y_spots
    if angle == 90:
        rotated["plot_x"] = y_spots - 1 - rotated["y"]
        rotated["plot_y"] = rotated["x"]
        return rotated, y_spots, x_spots
    if angle == 180:
        rotated["plot_x"] = x_spots - 1 - rotated["x"]
        rotated["plot_y"] = y_spots - 1 - rotated["y"]
        return rotated, x_spots, y_spots
    rotated["plot_x"] = rotated["y"]
    rotated["plot_y"] = x_spots - 1 - rotated["x"]
    return rotated, y_spots, x_spots


def summarize_lrs(frame: pd.DataFrame) -> pd.DataFrame:
    reads_frame = frame.copy()
    reads_frame["reads_num"] = pd.to_numeric(reads_frame["reads"], errors="coerce").fillna(0)
    summary = (
        reads_frame.groupby("LR", as_index=False)
        .agg(
            unique_SR=("SR", "nunique"),
            total_reads=("reads_num", "sum"),
            unique_UR=(
                "UR",
                lambda values: pd.Series(values)
                .astype(str)
                .str.strip()
                .replace("", pd.NA)
                .dropna()
                .nunique(),
            ),
        )
        .sort_values(
            ["unique_SR", "unique_UR", "total_reads", "LR"],
            ascending=[False, False, False, True],
            kind="stable",
        )
        .reset_index(drop=True)
    )
    return summary


def _plot_lr_spatial(
    cluster_frame: pd.DataFrame,
    lr_data: pd.DataFrame,
    output_path: Path,
    x_spots: int,
    y_spots: int,
    circle_color: str,
    rotate: int,
    cluster_alpha: float,
) -> None:
    cluster_frame, plot_x_spots, plot_y_spots = rotate_coordinates(
        cluster_frame,
        x_spots=x_spots,
        y_spots=y_spots,
        angle=rotate,
    )
    lr_spots = (
        lr_data.assign(_ur_clean=lr_data["UR"].astype(str).str.strip())
        .groupby(["SR", "LR", "x", "y"], as_index=False)
        .agg(
            unique_ur_count=(
                "_ur_clean",
                lambda values: values.replace("", pd.NA).dropna().nunique(),
            )
        )
        .sort_values(["SR", "x", "y"], kind="stable")
        .copy()
    )
    lr_spots, _, _ = rotate_coordinates(
        lr_spots,
        x_spots=x_spots,
        y_spots=y_spots,
        angle=rotate,
    )
    lr_spots["circle_size"] = lr_spots["unique_ur_count"].map(
        lambda count: 30 if count <= 1 else 70 if count == 2 else 140
    )
    mutations_text = "mutations unavailable"
    if "mutations" in lr_data.columns:
        mutations_values = (
            lr_data["mutations"]
            .astype(str)
            .str.strip()
            .replace("", pd.NA)
            .dropna()
            .drop_duplicates()
            .tolist()
        )
        if mutations_values:
            mutations_text = mutations_values[0]
    mutations_text = format_fixed_height_title(mutations_text)

    cluster_legend = (
        cluster_frame[["leiden", "color"]]
        .drop_duplicates()
        .sort_values("leiden", key=lambda col: col.astype(str))
        .reset_index(drop=True)
    )
    cluster_handles = [
        Patch(facecolor=row.color, edgecolor="none", alpha=cluster_alpha, label=f"{row.leiden}")
        for row in cluster_legend.itertuples(index=False)
    ]
    umi_bins = [(30, "1"), (70, "2"), (140, "3+")]
    umi_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="None",
            markerfacecolor="none",
            markeredgecolor=circle_color,
            markeredgewidth=2,
            markersize=legend_marker_size_from_scatter_area(area),
            label=label,
        )
        for area, label in umi_bins
    ]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, (ax, legend_ax) = plt.subplots(
        ncols=2,
        figsize=FIG_SIZE,
        gridspec_kw={"width_ratios": [5.7, 0.82]},
    )
    ax.scatter(
        cluster_frame["plot_x"],
        cluster_frame["plot_y"],
        s=28,
        c=cluster_frame["color"],
        marker="s",
        alpha=cluster_alpha,
        linewidths=0,
    )
    if not lr_spots.empty:
        ax.scatter(
            lr_spots["plot_x"],
            lr_spots["plot_y"],
            s=lr_spots["circle_size"],
            facecolors="none",
            edgecolors=circle_color,
            linewidths=1.5,
            clip_on=False,
        )

    ax.set_xlim(-PLOT_EDGE_PAD, plot_x_spots - 1 + PLOT_EDGE_PAD)
    ax.set_ylim(-PLOT_EDGE_PAD, plot_y_spots - 1 + PLOT_EDGE_PAD)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(mutations_text, fontsize=8, pad=5, loc="center")

    legend_ax.axis("off")
    cluster_legend_artist = legend_ax.legend(
        handles=cluster_handles,
        loc="upper left",
        frameon=False,
        fontsize=8,
        title="cluster",
        title_fontsize=9,
        borderaxespad=0,
        handletextpad=0.45,
        labelspacing=0.55,
        handlelength=1.15,
        handleheight=1.15,
        bbox_to_anchor=(0.0, 1.01),
    )
    legend_ax.add_artist(cluster_legend_artist)
    legend_ax.legend(
        handles=umi_handles,
        loc="upper left",
        frameon=False,
        fontsize=8,
        title="n_UMI",
        title_fontsize=9,
        borderaxespad=0,
        handletextpad=0.55,
        labelspacing=0.85,
        bbox_to_anchor=(0.0, 0.33),
    )

    fig.subplots_adjust(left=0.028, right=0.998, top=FIG_TOP, bottom=FIG_BOTTOM, wspace=0.0)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def build_file_stem(rank: int, unique_sr: int, unique_ur: int) -> str:
    return f"topLR_{rank:03d}_sr{unique_sr:03d}_ur{unique_ur:03d}"


def process_label(
    input_dir: Path,
    output_root: Path,
    label: str,
    cluster_frame: pd.DataFrame,
    top_n: int,
    rotate: int,
    cluster_alpha: float,
) -> None:
    input_file = resolve_input_path(input_dir, label)
    frame = pd.read_csv(input_file, dtype=str, keep_default_na=False)
    required_columns = {"SR", "UR", "LR", "reads", "x", "y"}
    missing_columns = required_columns - set(frame.columns)
    if missing_columns:
        raise SystemExit(
            f"Input file {input_file} is missing required columns: {missing_columns}"
        )

    frame = sanitize_table(frame)
    summary = summarize_lrs(frame).head(top_n).reset_index(drop=True)
    output_dir = output_root / label / "top_lr_plots"
    output_dir.mkdir(parents=True, exist_ok=True)

    x_spots = max(int(frame["x"].max()), int(cluster_frame["x"].max())) + 1
    y_spots = max(int(frame["y"].max()), int(cluster_frame["y"].max())) + 1
    circle_color = choose_contrasting_circle_color(cluster_frame)

    manifest_rows = []
    for idx, row in enumerate(summary.itertuples(index=False), start=1):
        lr = str(row.LR)
        lr_data = frame[frame["LR"] == lr].copy()
        stem = build_file_stem(idx, int(row.unique_SR), int(row.unique_UR))
        output_path = output_dir / f"{stem}.png"
        _plot_lr_spatial(
            cluster_frame=cluster_frame,
            lr_data=lr_data,
            output_path=output_path,
            x_spots=x_spots,
            y_spots=y_spots,
            circle_color=circle_color,
            rotate=rotate,
            cluster_alpha=cluster_alpha,
        )
        manifest_rows.append(
            {
                "rank": idx,
                "LR": lr,
                "unique_SR": int(row.unique_SR),
                "unique_UR": int(row.unique_UR),
                "total_reads": int(row.total_reads),
                "plot_file": output_path.name,
            }
        )
        print(f"Saved: {output_path}")

    manifest = pd.DataFrame(manifest_rows)
    manifest_path = output_dir / f"{label}_top_lr_plot_manifest.csv"
    manifest.to_csv(manifest_path, index=False)
    print(f"Saved: {manifest_path}")


def main() -> None:
    args = parse_args()
    input_dir = Path(args.input_dir)
    output_root = Path(args.output_dir) if args.output_dir else input_dir
    cluster_csv = Path(args.cluster_csv)
    cluster_frame = load_cluster_background(cluster_csv)
    output_root.mkdir(parents=True, exist_ok=True)

    for label in args.labels:
        print(f"\n=== Processing {label} ===")
        process_label(
            input_dir=input_dir,
            output_root=output_root,
            label=label,
            cluster_frame=cluster_frame,
            top_n=args.top_n,
            rotate=args.rotate,
            cluster_alpha=args.cluster_alpha,
        )


if __name__ == "__main__":
    main()
