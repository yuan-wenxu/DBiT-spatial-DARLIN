"""Shared I/O, spatial-coordinate, tissue-merge, and plotting helpers."""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from PIL import Image
from scipy.stats import pearsonr, spearmanr

DOWNSAMPLE_FACTOR = 10
MIN_OUTPUT_DIMENSION = 1500
PRIMARY_COLOR = "#0072B2"
HEATMAP_COLOR = "#D73027"
UMI_COLOR = HEATMAP_COLOR
GENE_COLOR = HEATMAP_COLOR
THRESHOLD_COLOR = "#D55E00"
CONTINUOUS_CMAP = "Reds"
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


def categorical_colors(count):
    return [
        CATEGORICAL_COLORS[index % len(CATEGORICAL_COLORS)]
        for index in range(count)
    ]


def alpha_colormap(color):
    rgba = np.ones((256, 4), dtype=float)
    rgba[:, :3] = mcolors.to_rgb(color)
    rgba[:, 3] = np.linspace(0, 1, 256)
    return mcolors.ListedColormap(rgba)


@dataclass(frozen=True)
class SpatialPlotConfig:
    x_spots_number: int
    y_spots_number: int
    length_spot: int
    interval: int
    pixel_length: float


@dataclass(frozen=True)
class ScatterConfig:
    xlabel: str
    ylabel: str
    title: str
    equal_axis: bool = False
    grid: bool = True
    saturation_x: float | None = None
    saturation_y: float | None = None


def open_text(path, mode="rt"):
    """Open plain-text or gzip-compressed input using text mode."""
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode, encoding="utf-8")


def read_nonempty_lines(path) -> list[str]:
    with open_text(path) as handle:
        return [line.strip() for line in handle if line.strip()]


def iter_fastq(handle):
    """Yield validated ``(read_id, sequence, quality)`` FASTQ records."""
    while True:
        id_line = handle.readline()
        if not id_line:
            break
        sequence_line = handle.readline()
        plus_line = handle.readline()
        quality_line = handle.readline()
        if not (sequence_line and plus_line and quality_line):
            raise ValueError("Incomplete FASTQ record encountered.")
        if not id_line.startswith("@") or not plus_line.startswith("+"):
            raise ValueError("Invalid FASTQ structure (missing @ or + line).")
        read_id = id_line[1:].strip()
        sequence = sequence_line.strip()
        quality = quality_line.strip()
        if len(sequence) != len(quality):
            raise ValueError(
                f"Length mismatch (seq {len(sequence)} vs qual {len(quality)}) "
                f"at read {read_id}"
            )
        yield read_id, sequence, quality


def iter_paired_fastq(handle1, handle2):
    """Yield validated records from two FASTQ files in lockstep."""
    iterator1 = iter_fastq(handle1)
    iterator2 = iter_fastq(handle2)
    while True:
        try:
            record1 = next(iterator1)
        except StopIteration:
            try:
                next(iterator2)
            except StopIteration:
                return
            raise ValueError("File 1 ended before file 2.")
        try:
            record2 = next(iterator2)
        except StopIteration as error:
            raise ValueError("File 2 ended before file 1.") from error
        yield (*record1, *record2)


def validate_barcode_length(cb_len: int) -> int:
    if (
        not isinstance(cb_len, int)
        or isinstance(cb_len, bool)
        or cb_len <= 0
        or cb_len % 2 != 0
    ):
        raise ValueError("cb_len must be a positive even integer")
    return cb_len // 2


def read_barcode_components(whitelist_path) -> list[str]:
    components = read_nonempty_lines(whitelist_path)
    if len(components) != len(set(components)):
        raise ValueError(f"Whitelist contains duplicate barcodes: {whitelist_path}")
    return components


def add_spatial_coordinates(
    data: pd.DataFrame,
    barcode_column: str,
    barcode_a_whitelist_path,
    cb_len: int,
    barcode_b_whitelist_path=None,
) -> pd.DataFrame:
    """Add x/y barcode components and zero-based grid coordinates."""
    component_len = validate_barcode_length(cb_len)
    barcode_as = read_barcode_components(barcode_a_whitelist_path)
    barcode_bs = read_barcode_components(
        barcode_b_whitelist_path or barcode_a_whitelist_path
    )
    x_coordinate = {barcode: index for index, barcode in enumerate(barcode_as)}
    y_coordinate = {barcode: index for index, barcode in enumerate(barcode_bs)}
    result = data.copy()
    barcodes = result[barcode_column].astype(str)
    result["xbc"] = barcodes.str[component_len:cb_len]
    result["ybc"] = barcodes.str[:component_len]
    result["x"] = result["xbc"].map(x_coordinate).fillna(-1).astype(int)
    result["y"] = result["ybc"].map(y_coordinate).fillna(-1).astype(int)
    return result


def load_cell_numbers(cell_number_file) -> pd.DataFrame:
    cells = pd.read_csv(cell_number_file)
    required = {"x", "y", "count", "in_tissue"}
    missing = required.difference(cells.columns)
    if missing:
        raise ValueError(f"Cell-count file is missing columns: {sorted(missing)}")
    if cells["in_tissue"].dtype != bool:
        cells["in_tissue"] = (
            cells["in_tissue"].astype(str).str.lower().isin(["true", "1", "yes"])
        )
    return cells


def merge_tissue_spots(data: pd.DataFrame, cells: pd.DataFrame) -> pd.DataFrame:
    merged = data.merge(
        cells[["x", "y", "count", "in_tissue"]], on=["x", "y"]
    )
    return merged[merged["in_tissue"]].copy()


def plot_spatial_heatmaps(
    table_path,
    barcode_a_whitelist_path,
    barcode_b_whitelist_path,
    output_path,
    cb_len: int,
    x_spots_number: int,
    y_spots_number: int,
) -> None:
    """Plot reads/UMI/gene metrics using separate A/B barcode whitelists."""
    data = pd.read_csv(table_path)
    if "x" not in data.columns or "y" not in data.columns:
        data = add_spatial_coordinates(
            data,
            "SR",
            barcode_a_whitelist_path,
            cb_len,
            barcode_b_whitelist_path=barcode_b_whitelist_path,
        )
    if "umi_count" not in data.columns and "UR" in data.columns:
        aggregations = {"UR": "nunique"}
        if "reads" in data.columns:
            aggregations["reads"] = "sum"
        data = data.groupby(["x", "y"]).agg(aggregations).reset_index()
        data["umi_count"] = data["UR"]

    full_index = pd.MultiIndex.from_product(
        [range(x_spots_number), range(y_spots_number)], names=["x", "y"]
    )
    metric_specs = (
        ("reads", "Reads counts", "Reads_counts_heatmap.png"),
        ("umi_count", "UMI counts", "UMI_counts_heatmap.png"),
        ("gene_count", "Gene counts", "Gene_counts_heatmap.png"),
    )
    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)
    for column, title, filename in metric_specs:
        if column not in data.columns:
            continue
        values = data[["x", "y", column]].copy()
        values[column] = pd.to_numeric(values[column], errors="coerce").fillna(0)
        values = values.groupby(["x", "y"])[column].sum().reindex(
            full_index, fill_value=0
        )
        pivot = values.unstack("x")
        figure, axis = plt.subplots(figsize=(5, 4))
        sns.heatmap(
            pivot,
            cmap=CONTINUOUS_CMAP,
            annot=False,
            linewidths=0,
            xticklabels=False,
            yticklabels=False,
            ax=axis,
        )
        axis.set(xlabel="", ylabel="", title=title)
        axis.title.set_fontsize(12)
        axis.xaxis.label.set_size(10)
        axis.yaxis.label.set_size(10)
        axis.tick_params(axis="both", labelsize=10)
        figure.tight_layout()
        figure.savefig(output_path / filename, dpi=300, bbox_inches="tight")
        plt.close(figure)


def plot_scatter(x, y, output_path, config: ScatterConfig) -> None:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if config.saturation_x:
        x = x / config.saturation_x
    if config.saturation_y:
        y = y / config.saturation_y

    figure, axis = plt.subplots(figsize=(5, 4))
    axis.scatter(x, y, s=10, alpha=0.1, color=PRIMARY_COLOR)
    if len(x) >= 2:
        pearson = pearsonr(x, y)
        spearman = spearmanr(x, y)
        axis.text(
            0.05,
            0.95,
            f"Pearson Correlation: {pearson.statistic:.3f} (p={pearson.pvalue:.3e})",
            transform=axis.transAxes,
            fontsize=8,
            va="top",
        )
        axis.text(
            0.05,
            0.90,
            f"Spearman Correlation: {spearman.statistic:.3f} (p={spearman.pvalue:.3e})",
            transform=axis.transAxes,
            fontsize=8,
            va="top",
        )
    if config.equal_axis:
        maximum = max(axis.get_xlim()[1], axis.get_ylim()[1])
        axis.set(xlim=(0, maximum), ylim=(0, maximum))
        axis.plot(
            [0, maximum],
            [0, maximum],
            linestyle="--",
            color=THRESHOLD_COLOR,
        )
    if config.grid:
        axis.grid(True)
    axis.set(xlabel=config.xlabel, ylabel=config.ylabel, title=config.title)
    axis.title.set_fontsize(12)
    axis.xaxis.label.set_size(10)
    axis.yaxis.label.set_size(10)
    axis.tick_params(axis="both", labelsize=10)
    figure.tight_layout()
    figure.savefig(Path(output_path) / f"{config.title}_scatter.png", dpi=300, bbox_inches="tight")
    plt.close(figure)


def spatial_image_size(config):
    width = int(
        (
            config.x_spots_number * config.length_spot
            + (config.x_spots_number - 1) * config.interval
        )
        / config.pixel_length
    )
    height = int(
        (
            config.y_spots_number * config.length_spot
            + (config.y_spots_number - 1) * config.interval
        )
        / config.pixel_length
    )
    return width, height


def resized_frame_size(config):
    original_width, original_height = spatial_image_size(config)
    return (
        max(MIN_OUTPUT_DIMENSION, original_width // DOWNSAMPLE_FACTOR),
        max(MIN_OUTPUT_DIMENSION, original_height // DOWNSAMPLE_FACTOR),
    )


def spatial_frame_shape(config):
    width = int(
        config.x_spots_number * config.length_spot
        + (config.x_spots_number - 1) * config.interval
    )
    height = int(
        config.y_spots_number * config.length_spot
        + (config.y_spots_number - 1) * config.interval
    )
    return height, width, 4


def scale_frame_alpha(frame):
    positive = frame[:, :, 3][frame[:, :, 3] > 0]
    if positive.size:
        p95 = np.percentile(positive, 95)
        frame[:, :, 3] = np.clip(frame[:, :, 3] / p95, 0, 1) * 255
    return frame.astype(np.uint8)


def plot_value_legend(data, output_path, name, color):
    values = np.asarray(data, dtype=float)
    values = values[np.isfinite(values) & (values > 0)]
    if values.size == 0:
        return
    p95 = np.percentile(values, 95)

    norm = mcolors.Normalize(vmin=0, vmax=p95)
    sm = plt.cm.ScalarMappable(cmap=alpha_colormap(color), norm=norm)
    sm.set_array([])
    fig = plt.figure(figsize=(5, 10))
    ax = fig.add_axes([0.10, 0.14, 0.14, 0.70])
    cbar = fig.colorbar(sm, cax=ax, orientation='vertical')
    ticks = np.linspace(0, p95, 5)
    cbar.set_ticks(ticks)
    if p95 >= 10:
        ticklabels = [f"{t:,.0f}" for t in ticks]
        p95_label = f"{p95:,.0f}"
    else:
        ticklabels = [f"{t:.2f}" for t in ticks]
        p95_label = f"{p95:.2f}"
    cbar.set_ticklabels(ticklabels, fontsize=22)
    cbar.ax.tick_params(labelsize=22, pad=10)
    display_name = name.replace('_', ' ')
    cbar.set_label(display_name, fontsize=22, labelpad=28)
    fig.suptitle(
        f"Color capped at P95\nP95 = {p95_label}",
        fontsize=22,
        y=0.96,
    )
    fig.text(
        0.5,
        0.04,
        "Values ≥ P95 use maximum opacity",
        ha='center',
        va='center',
        fontsize=16,
    )
    ax.tick_params(labelsize=22)
    fig.savefig(f'{output_path}/{name}', dpi=300)
    plt.close()

def plot_spatial_frames(data, output_path, config: SpatialPlotConfig):
    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)
    output_size = resized_frame_size(config)
    frame_shape = spatial_frame_shape(config)

    positive_cell_spots = data['count'] > 0 if 'count' in data.columns else None
    has_umi_per_cell = (
        'umi_count' in data.columns
        and positive_cell_spots is not None
        and bool((positive_cell_spots & (data['umi_count'] > 0)).any())
    )

    if 'umi_count' in data.columns:
        plot_value_legend(
            data['umi_count'], output_path, 'UMI_distribution', UMI_COLOR
        )
        if has_umi_per_cell:
            plot_value_legend(
                data.loc[positive_cell_spots, 'umi_count']
                / data.loc[positive_cell_spots, 'count'],
                output_path,
                'UMI_per_cell_distribution',
                UMI_COLOR,
            )
    if 'gene_count' in data.columns:
        plot_value_legend(
            data['gene_count'], output_path, 'Gene_distribution', GENE_COLOR
        )

    if 'leiden' in data.columns:
        frame = np.zeros(frame_shape, dtype=np.uint8)
        if 'color' in data.columns:
            cluster_colors = {
                int(row['leiden']): np.array(mcolors.to_rgba(row['color']))
                for _, row in data[['leiden', 'color']].drop_duplicates().iterrows()
            }
        else:
            clusters = sorted(data['leiden'].astype(int).unique())
            colors = categorical_colors(len(clusters))
            cluster_colors = {
                cluster_id: np.array(mcolors.to_rgba(colors[i]))
                for i, cluster_id in enumerate(clusters)
            }

    if 'umi_count' in data.columns:
        frame_umi = np.zeros(frame_shape, dtype=np.float32)
        umi_rgb = np.asarray(mcolors.to_rgb(UMI_COLOR)) * 255
        if has_umi_per_cell:
            frame_umi_per_cell = np.zeros(frame_shape, dtype=np.float32)

    if 'gene_count' in data.columns:
        frame_gene = np.zeros(frame_shape, dtype=np.float32)
        gene_rgb = np.asarray(mcolors.to_rgb(GENE_COLOR)) * 255

    for _, row in data.iterrows():
        x_idx = int(row['x'])
        y_idx = int(row['y'])

        x_start = x_idx * (config.length_spot + config.interval)
        y_start = y_idx * (config.length_spot + config.interval)
        x_end = x_start + config.length_spot
        y_end = y_start + config.length_spot

        if 'umi_count' in data.columns:
            frame_umi[y_start: y_end, x_start: x_end, :3] = umi_rgb
            frame_umi[y_start: y_end, x_start: x_end, 3] = row['umi_count']
            if has_umi_per_cell and row['count'] > 0:
                frame_umi_per_cell[
                    y_start: y_end, x_start: x_end, :3
                ] = umi_rgb
                frame_umi_per_cell[y_start: y_end, x_start: x_end, 3] = row['umi_count'] / row['count']
        if 'gene_count' in data.columns:
            frame_gene[y_start: y_end, x_start: x_end, :3] = gene_rgb
            frame_gene[y_start: y_end, x_start: x_end, 3] = int(row['gene_count'])
        if 'leiden' in data.columns:
            cluster_id = int(row['leiden'])
            frame[y_start: y_end, x_start: x_end, :] = (cluster_colors[cluster_id] * 255).astype(np.uint8)

    if 'leiden' in data.columns:
        img_umap = Image.fromarray(frame, mode = 'RGBA')
        img_umap = img_umap.resize(output_size, resample = Image.NEAREST)
        img_umap.save(output_path / 'umap_filtered.png')
    if 'umi_count' in data.columns:
        frame_umi = scale_frame_alpha(frame_umi)

        if has_umi_per_cell:
            frame_umi_per_cell = scale_frame_alpha(frame_umi_per_cell)

        img_umi = Image.fromarray(frame_umi, mode = 'RGBA')
        img_umi = img_umi.resize(output_size, resample = Image.NEAREST)
        img_umi.save(output_path / 'umi_filtered.png')

        if has_umi_per_cell:
            img_umi_per_cell = Image.fromarray(frame_umi_per_cell, mode = 'RGBA')
            img_umi_per_cell = img_umi_per_cell.resize(
                output_size, resample = Image.NEAREST
            )
            img_umi_per_cell.save(output_path / 'umi_per_cell_filtered.png')

    if 'gene_count' in data.columns:
        frame_gene = scale_frame_alpha(frame_gene)

        img_gene = Image.fromarray(frame_gene, mode = 'RGBA')
        img_gene = img_gene.resize(output_size, resample = Image.NEAREST)
        img_gene.save(output_path / 'gene_filtered.png')
