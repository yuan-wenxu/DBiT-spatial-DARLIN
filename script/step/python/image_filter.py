import argparse
import gzip
from dataclasses import dataclass
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from scipy.io import mmread, mmwrite


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter mRNA and DARLIN results to spots inside tissue.",
        add_help=False,
    )
    parser.add_argument(
        "--help",
        action="help",
        help="Show this help message and exit",
    )
    parser.add_argument(
        "--tissue_positions_file",
        required=True,
        type=Path,
        help="Compressed TSV containing tissue spot positions",
    )
    parser.add_argument(
        "--mrna_path",
        type=Path,
        help="mRNA GeneFull directory containing raw/spatial_metrics.csv",
    )
    parser.add_argument(
        "--darlin_path",
        type=Path,
        help="Directory containing CA, RA, and TA DARLIN results",
    )
    parser.add_argument(
        "--x_spots_number",
        type=int,
        default=50,
        help="Number of spot rows",
    )
    parser.add_argument(
        "--y_spots_number",
        type=int,
        default=50,
        help="Number of spot columns",
    )
    parser.add_argument(
        "--length_spot",
        type=int,
        default=20,
        help="Length of each spot in pixels",
    )
    parser.add_argument(
        "--interval",
        type=int,
        default=20,
        help="Interval between spots in pixels",
    )
    parser.add_argument(
        "--pixel_length",
        type=float,
        default=0.294,
        help="Length of each pixel in microns",
    )
    args = parser.parse_args()
    if args.mrna_path is None and args.darlin_path is None:
        parser.error("at least one of --mrna_path or --darlin_path is required")
    return args


@dataclass(frozen=True)
class SpatialPlotConfig:
    x_spots_number: int
    y_spots_number: int
    length_spot: int
    interval: int
    pixel_length: float


DARLIN_LOCI = ("CA", "RA", "TA")
DOWNSAMPLE_FACTOR = 10
MIN_OUTPUT_DIMENSION = 1500
SPATIAL_COLOR = "#D73027"
CATEGORICAL_COLORS = (
    "#0072B2",
    "#E69F00",
    "#009E73",
    "#CC79A7",
    "#56B4E9",
    "#D55E00",
    "#F0E442",
    "#222222",
    "#6A3D9A",
    "#1B9E77",
    "#E7298A",
    "#66A61E",
    "#E6AB02",
    "#A6761D",
    "#7570B3",
    "#A6CEE3",
    "#FDBF6F",
    "#B2DF8A",
    "#FB9A99",
    "#CAB2D6",
)


def load_tissue_positions(tissue_positions_file: Path) -> pd.DataFrame:
    positions = pd.read_csv(tissue_positions_file, sep="\t", compression="infer")
    required = {"array_row", "array_col"}
    missing = required.difference(positions.columns)
    if missing:
        raise ValueError(
            f"Tissue-position file is missing columns: {sorted(missing)}"
        )
    return positions.rename(columns={"array_row": "row", "array_col": "col"})


def merge_tissue_spots(
    data: pd.DataFrame, tissue_positions: pd.DataFrame
) -> pd.DataFrame:
    missing = {"row", "col"}.difference(data.columns)
    if missing:
        raise ValueError(f"Spatial table is missing columns: {sorted(missing)}")
    position_columns = [
        column
        for column in tissue_positions.columns
        if column not in {"barcode", "row", "col"}
    ]
    return data.merge(
        tissue_positions[["row", "col", *position_columns]],
        on=["row", "col"],
        how="inner",
    )


def mrna_matrix_output_path(mrna_path: Path) -> Path:
    for parent in mrna_path.parents:
        if parent.name == "results":
            return parent.parent / "matrix"
    raise ValueError(
        f"Cannot locate the mRNA directory from GeneFull path: {mrna_path}"
    )


def write_tissue_filtered_mrna_matrix(
    tissue_positions: pd.DataFrame,
    mrna_path: Path,
) -> Path:
    if "barcode" not in tissue_positions.columns:
        raise ValueError("Tissue-position file is missing column: barcode")

    raw_path = mrna_path / "raw"
    matrix_file = raw_path / "matrix.mtx"
    barcodes_file = raw_path / "barcodes.tsv"
    features_file = raw_path / "features.tsv"
    missing_files = [
        path.name
        for path in (matrix_file, barcodes_file, features_file)
        if not path.is_file()
    ]
    if missing_files:
        raise FileNotFoundError(
            f"Uncompressed 10x files missing from {raw_path}: {missing_files}"
        )

    matrix = mmread(matrix_file).tocsc()
    barcodes = pd.read_csv(
        barcodes_file,
        sep="\t",
        header=None,
        dtype=str,
        keep_default_na=False,
    )
    features = pd.read_csv(
        features_file,
        sep="\t",
        header=None,
        dtype=str,
        keep_default_na=False,
    )
    if matrix.shape != (len(features), len(barcodes)):
        raise ValueError(
            f"10x matrix shape {matrix.shape} does not match "
            f"{len(features)} features and {len(barcodes)} barcodes"
        )
    if barcodes.iloc[:, 0].duplicated().any():
        raise ValueError(f"Duplicate barcodes found in {barcodes_file}")

    tissue_barcodes = set(tissue_positions["barcode"].astype(str))
    retained_columns = [
        index
        for index, barcode in enumerate(barcodes.iloc[:, 0])
        if barcode in tissue_barcodes
    ]
    if not retained_columns:
        raise ValueError("No tissue-position barcodes matched the mRNA matrix")

    filtered_matrix = matrix[:, retained_columns].tocoo()
    filtered_barcodes = barcodes.iloc[retained_columns].reset_index(drop=True)
    output_path = mrna_matrix_output_path(mrna_path)
    output_path.mkdir(parents=True, exist_ok=True)
    with gzip.open(output_path / "matrix.mtx.gz", "wb") as handle:
        mmwrite(handle, filtered_matrix, field="integer")
    filtered_barcodes.to_csv(
        output_path / "barcodes.tsv.gz",
        sep="\t",
        header=False,
        index=False,
        compression="gzip",
    )
    features.to_csv(
        output_path / "features.tsv.gz",
        sep="\t",
        header=False,
        index=False,
        compression="gzip",
    )
    print(
        f"Tissue-filtered matrix spots: {filtered_matrix.shape[1]}/"
        f"{matrix.shape[1]}"
    )
    print(f"Wrote tissue-filtered 10x matrix: {output_path.resolve()}")
    return output_path


def categorical_colors(count: int) -> list[str]:
    return [
        CATEGORICAL_COLORS[index % len(CATEGORICAL_COLORS)]
        for index in range(count)
    ]


def alpha_colormap(color: str) -> mcolors.ListedColormap:
    rgba = np.ones((256, 4), dtype=float)
    rgba[:, :3] = mcolors.to_rgb(color)
    rgba[:, 3] = np.linspace(0, 1, 256)
    return mcolors.ListedColormap(rgba)


def spatial_frame_geometry(
    config: SpatialPlotConfig,
) -> tuple[tuple[int, int, int], tuple[int, int]]:
    width = (
        config.y_spots_number * config.length_spot
        + (config.y_spots_number - 1) * config.interval
    )
    height = (
        config.x_spots_number * config.length_spot
        + (config.x_spots_number - 1) * config.interval
    )
    full_width = int(width / config.pixel_length)
    full_height = int(height / config.pixel_length)
    output_size = (
        max(MIN_OUTPUT_DIMENSION, full_width // DOWNSAMPLE_FACTOR),
        max(MIN_OUTPUT_DIMENSION, full_height // DOWNSAMPLE_FACTOR),
    )
    return (height, width, 4), output_size


def scale_frame_alpha(frame: np.ndarray) -> np.ndarray:
    positive = frame[:, :, 3][frame[:, :, 3] > 0]
    if positive.size:
        percentile_95 = np.percentile(positive, 95)
        frame[:, :, 3] = np.clip(frame[:, :, 3] / percentile_95, 0, 1) * 255
    return frame.astype(np.uint8)


def plot_value_legend(
    data: pd.Series, output_path: Path, name: str, color: str
) -> None:
    values = np.asarray(data, dtype=float)
    values = values[np.isfinite(values) & (values > 0)]
    if values.size == 0:
        return
    percentile_95 = np.percentile(values, 95)
    norm = mcolors.Normalize(vmin=0, vmax=percentile_95)
    scalar_map = plt.cm.ScalarMappable(cmap=alpha_colormap(color), norm=norm)
    scalar_map.set_array([])
    figure = plt.figure(figsize=(5, 10))
    axis = figure.add_axes([0.10, 0.14, 0.14, 0.70])
    colorbar = figure.colorbar(scalar_map, cax=axis, orientation="vertical")
    ticks = np.linspace(0, percentile_95, 5)
    colorbar.set_ticks(ticks)
    if percentile_95 >= 10:
        tick_labels = [f"{tick:,.0f}" for tick in ticks]
        percentile_label = f"{percentile_95:,.0f}"
    else:
        tick_labels = [f"{tick:.2f}" for tick in ticks]
        percentile_label = f"{percentile_95:.2f}"
    colorbar.set_ticklabels(tick_labels, fontsize=22)
    colorbar.ax.tick_params(labelsize=22, pad=10)
    colorbar.set_label(name.replace("_", " "), fontsize=22, labelpad=28)
    figure.suptitle(
        f"Color capped at P95\nP95 = {percentile_label}", fontsize=22, y=0.96
    )
    figure.text(
        0.5,
        0.04,
        "Values ≥ P95 use maximum opacity",
        ha="center",
        va="center",
        fontsize=16,
    )
    figure.savefig(output_path / name, dpi=300)
    plt.close(figure)


def plot_spatial_frames(
    data: pd.DataFrame, output_path: Path, config: SpatialPlotConfig
) -> None:
    output_path.mkdir(parents=True, exist_ok=True)
    frame_shape, output_size = spatial_frame_geometry(config)
    spatial_rgb = np.asarray(mcolors.to_rgb(SPATIAL_COLOR)) * 255
    metric_frames = {}
    for column in ("umi_count", "gene_count"):
        if column not in data.columns:
            continue
        legend_name = {
            "umi_count": "UMI_distribution",
            "gene_count": "Gene_distribution",
        }[column]
        plot_value_legend(
            data[column], output_path, legend_name, SPATIAL_COLOR
        )
        frame = np.zeros(frame_shape, dtype=np.float32)
        metric_frames[column] = frame

    cluster_frame = None
    cluster_colors = {}
    if "leiden" in data.columns:
        cluster_frame = np.zeros(frame_shape, dtype=np.uint8)
        if "color" in data.columns:
            cluster_colors = {
                int(row["leiden"]): np.asarray(mcolors.to_rgba(row["color"]))
                for _, row in data[["leiden", "color"]].drop_duplicates().iterrows()
            }
        else:
            clusters = sorted(data["leiden"].astype(int).unique())
            colors = categorical_colors(len(clusters))
            cluster_colors = {
                cluster: np.asarray(mcolors.to_rgba(colors[index]))
                for index, cluster in enumerate(clusters)
            }

    for _, spot in data.iterrows():
        x_index = config.y_spots_number - 1 - int(spot["col"])
        y_index = int(spot["row"])
        x_start = x_index * (config.length_spot + config.interval)
        y_start = y_index * (config.length_spot + config.interval)
        x_end = x_start + config.length_spot
        y_end = y_start + config.length_spot
        for column, frame in metric_frames.items():
            frame[y_start:y_end, x_start:x_end, :3] = spatial_rgb
            frame[y_start:y_end, x_start:x_end, 3] = spot[column]
        if cluster_frame is not None:
            cluster = int(spot["leiden"])
            cluster_frame[y_start:y_end, x_start:x_end] = (
                cluster_colors[cluster] * 255
            ).astype(np.uint8)

    if cluster_frame is not None:
        image = Image.fromarray(cluster_frame, mode="RGBA")
        image.resize(output_size, resample=Image.Resampling.NEAREST).save(
            output_path / "umap_filtered.png"
        )
    for column, frame in metric_frames.items():
        image = Image.fromarray(scale_frame_alpha(frame), mode="RGBA")
        image.resize(output_size, resample=Image.Resampling.NEAREST).save(
            output_path / f"{column.split('_')[0]}_filtered.png"
        )


def filter_mrna_by_tissue(
    tissue_positions: pd.DataFrame,
    mrna_path: Path,
    frame_config: SpatialPlotConfig,
) -> Path:
    method_path = mrna_path / "raw"
    input_file = method_path / "spatial_metrics.csv"
    data = pd.read_csv(input_file)
    filtered = merge_tissue_spots(data, tissue_positions)
    output_file = method_path / "spatial_metrics_tissuefiltered.csv"
    filtered.to_csv(output_file, index=False)
    print(input_file)
    print(f"Tissue-filtered spots: {len(filtered)}")
    print(f'Total UMI: {filtered["umi_count"].sum()}')
    print(f'Total Gene: {filtered["gene_count"].sum()}')
    print(f'Mean UMI: {filtered["umi_count"].mean()}')
    print(f'Median UMI: {filtered["umi_count"].median()}')
    print(f'Mean Gene: {filtered["gene_count"].mean()}')
    print(f'Median Gene: {filtered["gene_count"].median()}')
    print()
    write_tissue_filtered_mrna_matrix(
        tissue_positions,
        mrna_path,
    )
    plot_spatial_frames(filtered, method_path, frame_config)
    return output_file


def filter_darlin_by_tissue(
    tissue_positions: pd.DataFrame, darlin_path: Path
) -> list[Path]:
    output_files = []
    for locus in DARLIN_LOCI:
        locus_path = darlin_path / locus
        input_file = locus_path / "final.csv"
        if not input_file.is_file():
            continue
        filtered = merge_tissue_spots(pd.read_csv(input_file), tissue_positions)
        output_file = locus_path / "tissuefiltered.csv"
        filtered.to_csv(output_file, index=False)
        output_files.append(output_file)
        print(f"{locus}: retained {len(filtered)} tissue-filtered rows")
        print(f"Wrote tissue-filtered DARLIN table: {output_file.resolve()}")
    if not output_files:
        raise FileNotFoundError(f"No CA, RA, or TA final.csv found in {darlin_path}")
    return output_files


def main():
    args = parse_args()
    tissue_positions = load_tissue_positions(args.tissue_positions_file)
    if args.mrna_path is not None:
        frame_config = SpatialPlotConfig(
            args.x_spots_number,
            args.y_spots_number,
            args.length_spot,
            args.interval,
            args.pixel_length,
        )
        filter_mrna_by_tissue(
            tissue_positions,
            args.mrna_path,
            frame_config,
        )
    if args.darlin_path is not None:
        filter_darlin_by_tissue(tissue_positions, args.darlin_path)


if __name__ == "__main__":
    main()
