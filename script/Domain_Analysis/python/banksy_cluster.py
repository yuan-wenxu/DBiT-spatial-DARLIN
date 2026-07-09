#!/usr/bin/env python3
"""Cluster RCTD cell-type weights into spatial domains with BANKSY."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import anndata as ad
from banksy.cluster_methods import run_Leiden_partition
from banksy.embed_banksy import generate_banksy_matrix
from banksy.initialize_banksy import initialize_banksy
from banksy_utils.umap_pca import pca_umap
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from PIL import Image, ImageOps


BASE_GRID_SPOTS = 50
BASE_PLOT_SIZE = 4.37
LEGEND_WIDTH = 0.63
VERTICAL_MARGIN = 0.73
SPOT_SIDE_LENGTH = 0.82
PLOT_EDGE_PAD = 0.9
DOWNSAMPLE_FACTOR = 10
MIN_OUTPUT_DIMENSION = 1500


def spatial_figure_layout(
    x_spots: int, y_spots: int
) -> tuple[tuple[float, float], list[float]]:
    """Scale the spatial panel beyond 50 spots while keeping the legend fixed."""
    inches_per_spot = BASE_PLOT_SIZE / BASE_GRID_SPOTS
    plot_width = max(BASE_PLOT_SIZE, x_spots * inches_per_spot)
    plot_height = max(BASE_PLOT_SIZE, y_spots * inches_per_spot)
    return (
        (plot_width + LEGEND_WIDTH, plot_height + VERTICAL_MARGIN),
        [plot_width, LEGEND_WIDTH],
    )


def draw_spot_grid(axis, frame: pd.DataFrame, colors) -> None:
    """Draw equal-sized grid cells with a uniform gap in data coordinates."""
    half_side = SPOT_SIDE_LENGTH / 2
    cells = [
        Rectangle(
            (row.plot_x - half_side, row.plot_y - half_side),
            SPOT_SIDE_LENGTH,
            SPOT_SIDE_LENGTH,
        )
        for row in frame.itertuples(index=False)
    ]
    collection = PatchCollection(
        cells,
        facecolors=list(colors),
        edgecolors="none",
    )
    axis.add_collection(collection)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Cluster the cell_type_weights.csv output from spacexr.R "
            "into spatial domains with BANKSY and Leiden."
        )
    )
    parser.add_argument(
        "-i",
        "--weights-file",
        required=True,
        type=Path,
        help="RCTD cell_type_weights.csv file.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        help=(
            "New or empty output directory (default: banksy_output next to "
            "the input file)."
        ),
    )
    parser.add_argument(
        "--lambda-param",
        type=float,
        default=0.8,
        help="Contribution of spatial-neighbour features (default: 0.8).",
    )
    parser.add_argument(
        "--resolution",
        type=float,
        default=1.0,
        help="Leiden clustering resolution (default: 1.0).",
    )
    parser.add_argument(
        "--spatial-neighbors",
        type=int,
        default=30,
        help="Number of spatial neighbours used by BANKSY (default: 30).",
    )
    parser.add_argument(
        "--cluster-neighbors",
        type=int,
        default=50,
        help="Number of neighbours used for Leiden clustering (default: 50).",
    )
    parser.add_argument(
        "--pca-components",
        type=int,
        default=20,
        help="Number of BANKSY principal components (default: 20).",
    )
    parser.add_argument(
        "--max-m",
        type=int,
        choices=(0, 1),
        default=1,
        help="Maximum azimuthal transform order (default: 1).",
    )
    parser.add_argument(
        "--neighbor-decay",
        choices=("scaled_gaussian", "reciprocal"),
        default="scaled_gaussian",
        help="Spatial-neighbour weight decay (default: scaled_gaussian).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Leiden partition seed (default: 42).",
    )
    parser.add_argument(
        "--x-spots-number",
        type=int,
        required=True,
        help="Number of spots along the x axis in the complete DBiT grid.",
    )
    parser.add_argument(
        "--y-spots-number",
        type=int,
        required=True,
        help="Number of spots along the y axis in the complete DBiT grid.",
    )
    parser.add_argument(
        "--length-spot",
        type=int,
        required=True,
        help="Spot square size in pixels before resizing (default: 20).",
    )
    parser.add_argument(
        "--interval",
        type=int,
        required=True,
        help="Gap between neighbouring spots in pixels before resizing (default: 20).",
    )
    parser.add_argument(
        "--pixel-length",
        type=float,
        required=True,
        help="Pixel scaling factor used by frame_filtered.py (default: 1.0).",
    )
    parser.add_argument(
        "--orientation",
        choices=("normal", "horizontal", "vertical", "rotate"),
        help="Grid origin orientation: normal, horizontal, vertical, or rotate.",
    )
    parser.add_argument(
        "--swap_xy",
        action="store_true",
        help="Swap x and y axes after applying orientation.",
    )
    parser.add_argument(
        "--rotate",
        type=int,
        choices=(0, 90, 180, 270),
        default=0,
        help="Rotate the plotted spatial grid clockwise (default: 0).",
    )
    parser.add_argument(
        "--no-normalize-weights",
        action="store_true",
        help="Do not normalize each spot's RCTD weights to sum to one.",
    )
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if not args.weights_file.is_file():
        raise SystemExit(f"Weights file not found: {args.weights_file}")
    if args.output_dir is None:
        args.output_dir = args.weights_file.parent / "banksy_output"
    if not 0.0 <= args.lambda_param <= 1.0:
        raise SystemExit("--lambda-param must be between 0 and 1")
    if args.resolution <= 0:
        raise SystemExit("--resolution must be greater than zero")
    if args.spatial_neighbors < 1 or args.cluster_neighbors < 1:
        raise SystemExit("Neighbour counts must be positive integers")
    if args.pca_components < 1:
        raise SystemExit("--pca-components must be a positive integer")
    if args.x_spots_number < 1 or args.y_spots_number < 1:
        raise SystemExit("x/y spots number must be positive integers")
    if args.length_spot < 1:
        raise SystemExit("--length-spot must be a positive integer")
    if args.interval < 0:
        raise SystemExit("--interval must be a non-negative integer")
    if args.pixel_length <= 0:
        raise SystemExit("--pixel-length must be greater than zero")
    args.output_dir.mkdir(parents=True, exist_ok=True)


def load_weights(path: Path, normalize: bool):
    table = pd.read_csv(path)
    metadata_columns = ["barcode", "x", "y"]
    missing = [column for column in metadata_columns if column not in table.columns]
    if missing:
        raise SystemExit(
            "Weights file is missing required column(s): " + ", ".join(missing)
        )

    feature_columns = [column for column in table.columns if column not in metadata_columns]
    if not feature_columns:
        raise SystemExit("Weights file does not contain any cell-type columns")
    if table["barcode"].isna().any() or table["barcode"].astype(str).duplicated().any():
        raise SystemExit("barcode values must be non-empty and unique")

    numeric_columns = ["x", "y", *feature_columns]
    try:
        numeric = table[numeric_columns].apply(pd.to_numeric, errors="raise")
    except (TypeError, ValueError) as error:
        raise SystemExit(f"Coordinates and cell-type weights must be numeric: {error}") from error
    values = numeric[feature_columns].to_numpy(dtype=np.float64)
    coordinates = numeric[["x", "y"]].to_numpy(dtype=np.float64)
    if not np.isfinite(coordinates).all():
        raise SystemExit("Spatial coordinates contain missing or non-finite values")
    if not np.isfinite(values).all() or (values < 0).any():
        raise SystemExit("Cell-type weights must be finite and non-negative")

    row_sums = values.sum(axis=1)
    if (row_sums <= 0).any():
        bad_barcodes = table.loc[row_sums <= 0, "barcode"].astype(str).tolist()
        preview = ", ".join(bad_barcodes[:5])
        raise SystemExit(f"Cell-type weights sum to zero for spot(s): {preview}")
    if normalize:
        values = values / row_sums[:, None]

    barcodes = table["barcode"].astype(str).to_numpy()
    obs = pd.DataFrame(
        {"x": coordinates[:, 0], "y": coordinates[:, 1]}, index=barcodes
    )
    obs.index.name = "barcode"
    var = pd.DataFrame(index=pd.Index(feature_columns, name="cell_type"))
    adata = ad.AnnData(X=values, obs=obs, var=var)
    adata.obsm["spatial"] = coordinates
    return adata, table, feature_columns


def transform_coordinates(
    frame,
    x_spots: int,
    y_spots: int,
    orientation: str,
    swap_xy: bool,
    rotate: int,
):
    transformed = frame.copy()
    if orientation == "horizontal":
        transformed["x"] = x_spots - 1 - transformed["x"]
    elif orientation == "vertical":
        transformed["y"] = y_spots - 1 - transformed["y"]
    elif orientation == "rotate":
        transformed["x"] = x_spots - 1 - transformed["x"]
        transformed["y"] = y_spots - 1 - transformed["y"]
    if swap_xy:
        transformed["x"], transformed["y"] = (
            transformed["y"].copy(),
            transformed["x"].copy(),
        )
        x_spots, y_spots = y_spots, x_spots
    if rotate == 90:
        old_x = transformed["x"].copy()
        transformed["x"] = y_spots - 1 - transformed["y"]
        transformed["y"] = old_x
        x_spots, y_spots = y_spots, x_spots
    elif rotate == 180:
        transformed["x"] = x_spots - 1 - transformed["x"]
        transformed["y"] = y_spots - 1 - transformed["y"]
    elif rotate == 270:
        old_x = transformed["x"].copy()
        transformed["x"] = transformed["y"]
        transformed["y"] = x_spots - 1 - old_x
        x_spots, y_spots = y_spots, x_spots
    transformed["plot_x"] = transformed["x"]
    transformed["plot_y"] = transformed["y"]
    return transformed, x_spots, y_spots


def spot_frame_size(
    x_spots: int,
    y_spots: int,
    length_spot: int,
    interval: int,
) -> tuple[int, int]:
    width = int(x_spots * length_spot + (x_spots - 1) * interval)
    height = int(y_spots * length_spot + (y_spots - 1) * interval)
    return width, height


def render_spot_frame(
    frame: pd.DataFrame,
    color_map: dict,
    x_spots: int,
    y_spots: int,
    length_spot: int,
    interval: int,
) -> np.ndarray:
    width, height = spot_frame_size(x_spots, y_spots, length_spot, interval)
    canvas = np.zeros((height, width, 4), dtype=np.uint8)
    for row in frame.itertuples(index=False):
        x_idx = int(row.x)
        y_idx = int(row.y)
        x_start = x_idx * (length_spot + interval)
        y_start = y_idx * (length_spot + interval)
        x_end = x_start + length_spot
        y_end = y_start + length_spot
        canvas[y_start:y_end, x_start:x_end, :] = color_map[row.cluster]
    return canvas


def resize_spot_image(
    image: Image.Image,
    x_spots: int,
    y_spots: int,
    length_spot: int,
    interval: int,
    pixel_length: float,
) -> Image.Image:
    width, height = spot_frame_size(x_spots, y_spots, length_spot, interval)
    resized_width = max(
        MIN_OUTPUT_DIMENSION,
        int(width / pixel_length) // DOWNSAMPLE_FACTOR,
    )
    resized_height = max(
        MIN_OUTPUT_DIMENSION,
        int(height / pixel_length) // DOWNSAMPLE_FACTOR,
    )
    return image.resize((resized_width, resized_height), resample=Image.NEAREST)


def run_banksy(adata, args: argparse.Namespace):
    if adata.n_obs < 3:
        raise SystemExit("BANKSY clustering requires at least three spots")

    spatial_neighbors = min(args.spatial_neighbors, adata.n_obs - 1)
    banksy_dict = initialize_banksy(
        adata,
        coord_keys=("x", "y", "spatial"),
        num_neighbours=spatial_neighbors,
        nbr_weight_decay=args.neighbor_decay,
        max_m=args.max_m,
        plt_edge_hist=False,
        plt_nbr_weights=False,
        plt_agf_angles=False,
        plt_theta=False,
    )
    banksy_dict, banksy_adata = generate_banksy_matrix(
        adata,
        banksy_dict,
        lambda_list=[float(args.lambda_param)],
        max_m=args.max_m,
        plot_std=False,
        save_matrix=False,
        verbose=False,
    )

    max_components = min(banksy_adata.n_obs - 1, banksy_adata.n_vars)
    pca_components = min(args.pca_components, max_components)
    pca_umap(
        banksy_dict,
        pca_dims=[pca_components],
        plt_remaining_var=False,
        add_umap=False,
    )
    cluster_neighbors = min(args.cluster_neighbors, adata.n_obs - 1)
    results, _ = run_Leiden_partition(
        banksy_dict=banksy_dict,
        resolutions=[float(args.resolution)],
        num_nn=cluster_neighbors,
        num_iterations=-1,
        partition_seed=args.seed,
        match_labels=False,
        annotations=None,
    )
    if len(results) != 1:
        raise RuntimeError(f"Expected one BANKSY result, received {len(results)}")

    result = results.iloc[0]
    labels = result["labels"]
    labels = labels.dense if hasattr(labels, "dense") else labels
    result_adata = result["adata"]
    result_adata.obsm["spatial"] = result_adata.obs[["x", "y"]].to_numpy()
    result_adata.obs["banksy_cluster"] = [str(label) for label in labels]
    result_adata.obs["banksy_cluster"] = result_adata.obs["banksy_cluster"].astype(
        "category"
    )
    return result_adata, labels, spatial_neighbors, cluster_neighbors, pca_components


def save_outputs(
    adata,
    source_table,
    feature_columns: list[str],
    labels,
    args: argparse.Namespace,
    spatial_neighbors: int,
    cluster_neighbors: int,
    pca_components: int,
) -> None:
    clusters = source_table.copy()
    clusters["banksy_cluster"] = labels
    clusters.to_csv(args.output_dir / "banksy_clusters.csv", index=False)
    adata.write_h5ad(args.output_dir / "banksy_result.h5ad")

    plot_table = pd.DataFrame(
        {
            "x": adata.obs["x"].to_numpy(),
            "y": adata.obs["y"].to_numpy(),
            "cluster": labels,
        }
    )
    coordinates = plot_table[["x", "y"]].to_numpy(dtype=float)
    if not np.allclose(coordinates, np.rint(coordinates)):
        raise SystemExit("Square grid output requires integer x/y coordinates")
    plot_table[["x", "y"]] = np.rint(coordinates).astype(int)
    outside_grid = (
        (plot_table["x"] < 0)
        | (plot_table["x"] >= args.x_spots_number)
        | (plot_table["y"] < 0)
        | (plot_table["y"] >= args.y_spots_number)
    )
    if outside_grid.any():
        invalid = plot_table.loc[outside_grid, ["x", "y"]].iloc[0]
        raise SystemExit(
            f"Coordinate ({invalid['x']}, {invalid['y']}) is outside the configured "
            f"{args.x_spots_number}x{args.y_spots_number} grid"
        )
    grid_table = plot_table.copy()
    plot_table, plot_x_spots, plot_y_spots = transform_coordinates(
        plot_table,
        x_spots=args.x_spots_number,
        y_spots=args.y_spots_number,
        orientation=args.orientation,
        swap_xy=args.swap_xy,
        rotate=args.rotate,
    )

    cluster_ids = sorted(pd.unique(plot_table["cluster"]), key=int)
    cmap = plt.get_cmap("tab20", len(cluster_ids))
    cluster_colors = {
        cluster_id: np.rint(np.asarray(cmap(index)) * 255).astype(np.uint8)
        for index, cluster_id in enumerate(cluster_ids)
    }
    legend_colors = {
        cluster_id: mcolors.to_hex(cluster_colors[cluster_id] / 255.0)
        for cluster_id in cluster_ids
    }
    point_colors = plot_table["cluster"].map(legend_colors)

    figure_size, width_ratios = spatial_figure_layout(plot_x_spots, plot_y_spots)
    figure, (axis, legend_axis) = plt.subplots(
        ncols=2,
        figsize=figure_size,
        gridspec_kw={"width_ratios": width_ratios},
    )
    draw_spot_grid(axis, plot_table, point_colors)
    axis.set_xlim(-PLOT_EDGE_PAD, plot_x_spots - 1 + PLOT_EDGE_PAD)
    axis.set_ylim(-PLOT_EDGE_PAD, plot_y_spots - 1 + PLOT_EDGE_PAD)
    axis.set_aspect("equal")
    axis.invert_yaxis()
    axis.set_xticks([])
    axis.set_yticks([])
    for spine in axis.spines.values():
        spine.set_visible(False)
    axis.set_title("BANKSY spatial domains")
    legend_handles = [
        plt.Line2D(
            [0],
            [0],
            marker="s",
            linestyle="None",
            markerfacecolor=legend_colors[cluster_id],
            markeredgewidth=0,
            label=str(cluster_id),
        )
        for cluster_id in cluster_ids
    ]
    legend_axis.axis("off")
    legend_axis.legend(
        handles=legend_handles,
        title="cluster",
        loc="upper left",
        frameon=False,
        fontsize=8,
        title_fontsize=9,
        borderaxespad=0,
        handletextpad=0.45,
        labelspacing=0.55,
        handlelength=1.15,
        handleheight=1.15,
        bbox_to_anchor=(0.0, 1.01),
    )
    figure.subplots_adjust(left=0.028, right=0.998, top=0.92, bottom=0.04, wspace=0.0)
    figure.savefig(args.output_dir / "banksy_clusters.png", dpi=300)
    plt.close(figure)

    grid_dir = args.output_dir / "cluster_grids"
    grid_dir.mkdir(parents=True, exist_ok=True)
    for cluster_id in cluster_ids:
        members = grid_table[grid_table["cluster"] == cluster_id]
        grid = render_spot_frame(
            members[["x", "y", "cluster"]],
            color_map=cluster_colors,
            x_spots=args.x_spots_number,
            y_spots=args.y_spots_number,
            length_spot=args.length_spot,
            interval=args.interval,
        )
        resize_spot_image(
            Image.fromarray(grid, mode="RGBA"),
            x_spots=args.x_spots_number,
            y_spots=args.y_spots_number,
            length_spot=args.length_spot,
            interval=args.interval,
            pixel_length=args.pixel_length,
        ).save(grid_dir / f"cluster_{int(cluster_id):03d}.png")

    parameters = {
        "weights_file": str(args.weights_file.resolve()),
        "normalized_weights": not args.no_normalize_weights,
        "cell_type_columns": feature_columns,
        "lambda_param": args.lambda_param,
        "resolution": args.resolution,
        "spatial_neighbors": spatial_neighbors,
        "cluster_neighbors": cluster_neighbors,
        "pca_components": pca_components,
        "max_m": args.max_m,
        "neighbor_decay": args.neighbor_decay,
        "seed": args.seed,
        "x_spots_number": args.x_spots_number,
        "y_spots_number": args.y_spots_number,
        "length_spot": args.length_spot,
        "interval": args.interval,
        "pixel_length": args.pixel_length,
        "orientation": args.orientation,
        "swap_xy": args.swap_xy,
        "rotate": args.rotate,
    }
    with (args.output_dir / "run_parameters.json").open("w", encoding="utf-8") as handle:
        json.dump(parameters, handle, indent=2, ensure_ascii=False)
        handle.write("\n")


def main() -> None:
    args = parse_args()
    validate_args(args)
    adata, source_table, feature_columns = load_weights(
        args.weights_file, normalize=not args.no_normalize_weights
    )
    result = run_banksy(adata, args)
    save_outputs(
        result[0],
        source_table,
        feature_columns,
        result[1],
        args,
        result[2],
        result[3],
        result[4],
    )
    print(f"BANKSY clustering completed: {args.output_dir}")


if __name__ == "__main__":
    main()
