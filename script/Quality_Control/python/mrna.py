import argparse
from pathlib import Path

import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex, to_rgba
from matplotlib.patches import Patch
import numpy as np
from PIL import Image, ImageDraw
from scipy import sparse
from sklearn.neighbors import NearestNeighbors
from scipy.io import mmread
import seaborn as sns

from utils import (
    SpatialPlotConfig,
    add_spatial_coordinates,
    plot_spatial_heatmaps,
    spatial_frame_shape,
    spatial_image_size,
    validate_barcode_length,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot mRNA data",
        add_help=False,
    )
    parser.add_argument(
        "--help",
        action="help",
        help="Show this help message and exit",
    )
    parser.add_argument(
        "--file_path",
        required=True,
        help="Path to directory containing GeneFull or Gene folder",
    )
    parser.add_argument(
        "--barcodeA_whitelist",
        required=True,
        help="TSV/text file with one barcode A sequence per line",
    )
    parser.add_argument(
        "--barcodeB_whitelist",
        required=True,
        help="TSV/text file with one barcode B sequence per line",
    )
    parser.add_argument(
        "--cb-len",
        type=int,
        required=True,
        help="Total concatenated cell-barcode length in bp",
    )
    parser.add_argument(
        "--umi_min",
        type=int,
        default=900,
        help="Minimum UMI count per spot",
    )
    parser.add_argument(
        "--gene_min",
        type=int,
        default=300,
        help="Minimum gene count per spot",
    )
    parser.add_argument(
        "--min_cells",
        type=int,
        default=3,
        help="Minimum number of cells per gene",
    )
    parser.add_argument(
        "--x_spots_number",
        type=int,
        default=50,
        help="Number of spots in x direction",
    )
    parser.add_argument(
        "--y_spots_number",
        type=int,
        default=50,
        help="Number of spots in y direction",
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
    return parser.parse_args()


RANDOM_STATE = 42
AXIS_FONT_SIZE = 10
TITLE_FONT_SIZE = 12
UMI_COLOR = "#0072B2"
GENE_COLOR = "#E69F00"
THRESHOLD_COLOR = "#D55E00"
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


def set_axis_font_sizes(axis):
    """Apply the common title, axis-label, and tick font sizes."""
    axis.title.set_fontsize(TITLE_FONT_SIZE)
    axis.xaxis.label.set_size(AXIS_FONT_SIZE)
    axis.yaxis.label.set_size(AXIS_FONT_SIZE)
    axis.tick_params(axis="both", labelsize=AXIS_FONT_SIZE)


def load_and_filter_counts(
    file_path, apply_filter, umi_min=900, gene_min=300, min_cells=3
):
    """Load a STARsolo matrix, write QC plots, and apply spot/gene filters."""
    file_path = Path(file_path)
    matrix = mmread(file_path / "matrix.mtx")
    barcodes = pd.read_csv(
        file_path / "barcodes.tsv", header=None, sep="\t", dtype=str
    )
    features = pd.read_csv(
        file_path / "features.tsv", header=None, sep="\t", dtype=str
    )
    adata = sc.AnnData(X=matrix.T.tocsr())
    adata.obs_names = barcodes.iloc[:, 0].astype(str).to_numpy()
    adata.var_names = features.iloc[:, 0].astype(str).to_numpy()
    if features.shape[1] > 1:
        adata.var["gene_name"] = pd.array(features.iloc[:, 1], dtype="string")
    if features.shape[1] > 2:
        adata.var["feature_type"] = features.iloc[:, 2].astype(str).to_numpy()
    adata.obs_names.name = "barcode"
    adata.var_names.name = "gene"
    sc.pp.calculate_qc_metrics(adata, inplace=True)

    print(f"Before filtering: {adata.n_obs} spots, {adata.n_vars} genes")
    if not apply_filter:
        return adata

    plot_count_histogram(
        adata.obs["n_genes_by_counts"].to_numpy(),
        gene_min,
        "Gene Minimum",
        "Genes per Spot",
        "Number of Genes",
        file_path / "gene_counts_hist.png",
        color=GENE_COLOR,
    )
    plot_count_histogram(
        adata.obs["total_counts"].to_numpy(),
        umi_min,
        "UMI Minimum",
        "UMI per Spot",
        "Number of UMIs",
        file_path / "umi_counts_hist.png",
        color=UMI_COLOR,
    )
    spots_per_gene = np.asarray((adata.X > 0).sum(axis=0)).ravel()
    plot_spots_per_gene(spots_per_gene, adata.n_obs, min_cells, file_path)

    adata = adata[
        (adata.obs["total_counts"] >= umi_min)
        & (adata.obs["n_genes_by_counts"] >= gene_min),
        :,
    ].copy()
    sc.pp.filter_genes(adata, min_cells=min_cells)
    print(f"After filtering: {adata.n_obs} spots, {adata.n_vars} genes")
    plot_count_violins(adata, file_path / "violin_filtered_.png")
    return adata


def plot_count_violins(adata, output_path):
    figure, axes = plt.subplots(1, 2, figsize=(10, 4))
    sns.violinplot(
        y=adata.obs["n_genes_by_counts"],
        ax=axes[0],
        color=GENE_COLOR,
        inner="box",
        width=0.8,
    )
    axes[0].set(title="Genes per Spot", ylabel="Number of Genes")
    sns.violinplot(
        y=adata.obs["total_counts"],
        ax=axes[1],
        color=UMI_COLOR,
        inner="box",
        width=0.8,
    )
    axes[1].set(title="UMI per Spot", ylabel="Number of UMIs")
    for axis in axes:
        set_axis_font_sizes(axis)
    figure.tight_layout()
    figure.savefig(output_path, bbox_inches="tight", dpi=300)
    plt.close(figure)


def plot_count_histogram(
    values, cutoff, cutoff_label, title, xlabel, output_path, color
):
    p95 = np.percentile(values, 95)
    figure, axis = plt.subplots(figsize=(5, 4))
    axis.hist(values[values <= p95], bins=100, color=color)
    axis.axvline(
        cutoff,
        color=THRESHOLD_COLOR,
        linestyle="--",
        label=f"{cutoff_label}: {cutoff}",
    )
    axis.set_xlim(0, p95)
    axis.set_title(
        f"{title}\nUpper 5% omitted (P95 = {p95:,.0f})",
    )
    axis.set_xlabel(xlabel)
    axis.set_ylabel("Frequency")
    set_axis_font_sizes(axis)
    axis.legend(loc="upper right")
    figure.tight_layout()
    figure.savefig(output_path, bbox_inches="tight", dpi=300)
    plt.close(figure)


def plot_spots_per_gene(values, n_spots, min_cells, output_dir):
    bin_width = max(1, int(np.ceil(n_spots / 100)))
    plot_specs = (
        (
            np.arange(0, n_spots + bin_width, bin_width),
            None,
            output_dir / "spots_per_gene_hist.png",
        ),
        (
            np.arange(0.5, 21.5, 1),
            (0.5, 20.5),
            output_dir / "spots_per_gene_hist_small.png",
        ),
    )
    for bins, xlim, output_path in plot_specs:
        figure, axis = plt.subplots(figsize=(5, 4))
        axis.hist(values, bins=bins, color=GENE_COLOR)
        axis.axvline(
            min_cells,
            color=THRESHOLD_COLOR,
            linestyle="--",
            label=f"Minimum Cells: {min_cells}",
        )
        axis.set(title="Spots per Gene", xlabel="Number of Spots", ylabel="Frequency")
        if xlim:
            axis.set_xlim(*xlim)
            axis.set_xticks(range(1, 21))
        set_axis_font_sizes(axis)
        axis.legend(loc="upper right")
        figure.tight_layout()
        figure.savefig(output_path, bbox_inches="tight", dpi=300)
        plt.close(figure)


def sct_normalize_and_pca(adata, n_top_genes=3000, n_comps=50, theta=100):
    n_top_genes = min(n_top_genes, adata.n_vars)

    # Low-depth spots can have counts only in genes excluded by HVG selection.
    # Remove those zero-sum HVG rows before Pearson residual normalization,
    # which would otherwise produce NaN values from 0/0.
    while True:
        sc.experimental.pp.highly_variable_genes(
            adata,
            flavor="pearson_residuals",
            n_top_genes=n_top_genes,
            theta=theta,
            inplace=True,
        )
        hvg_mask = adata.var["highly_variable"].to_numpy()
        hvg_totals = np.asarray(adata[:, hvg_mask].X.sum(axis=1)).ravel()
        keep_spots = hvg_totals > 0
        removed_spots = int((~keep_spots).sum())
        if removed_spots == 0:
            break
        print(
            f"Removing {removed_spots} spots with zero counts across "
            "selected highly variable genes."
        )
        adata._inplace_subset_obs(keep_spots)

    n_comps = min(n_comps, adata.n_obs - 1, adata.n_vars - 1)
    if n_comps < 1:
        raise ValueError("Not enough observations or genes for PCA.")

    sc.experimental.pp.recipe_pearson_residuals(
        adata,
        theta=theta,
        n_top_genes=n_top_genes,
        n_comps=n_comps,
        random_state=RANDOM_STATE,
        inplace=True,
    )
    return n_comps


def build_snn_graph(adata, n_neighbors=30, n_pcs=20):
    if "X_pca" not in adata.obsm:
        raise ValueError("PCA coordinates not found in adata.obsm['X_pca'].")

    n_obs = adata.n_obs
    n_neighbors = min(n_neighbors, n_obs - 1)
    n_pcs = min(n_pcs, adata.obsm["X_pca"].shape[1])
    if n_neighbors < 1:
        raise ValueError("Not enough observations to construct an SNN graph.")

    x_pca = adata.obsm["X_pca"][:, :n_pcs]
    knn = NearestNeighbors(n_neighbors=n_neighbors + 1, metric="euclidean")
    knn.fit(x_pca)
    distances, indices = knn.kneighbors(x_pca)
    neighbor_indices = indices[:, 1:]
    neighbor_distances = distances[:, 1:]

    rows = np.repeat(np.arange(n_obs), n_neighbors)
    cols = neighbor_indices.ravel()
    data = np.ones(rows.shape[0], dtype=np.float32)
    knn_graph = sparse.csr_matrix((data, (rows, cols)), shape=(n_obs, n_obs))

    shared = knn_graph @ knn_graph.T
    shared.setdiag(0)
    shared.eliminate_zeros()
    shared = shared.tocoo()
    weights = shared.data / (2 * n_neighbors - shared.data)
    connectivities = sparse.csr_matrix((weights, (shared.row, shared.col)), shape=(n_obs, n_obs))
    connectivities = connectivities.maximum(connectivities.T)

    distances_graph = sparse.csr_matrix((neighbor_distances.ravel(), (rows, cols)), shape=(n_obs, n_obs))
    distances_graph = distances_graph.maximum(distances_graph.T)

    adata.obsp["distances"] = distances_graph
    adata.obsp["connectivities"] = connectivities
    adata.uns["neighbors"] = {
        "connectivities_key": "connectivities",
        "distances_key": "distances",
        "params": {
            "n_neighbors": n_neighbors,
            "n_pcs": n_pcs,
            "method": "umap",
            "graph_type": "snn",
            "snn_weight": "jaccard",
            "metric": "euclidean",
        },
    }


def make_h5ad_names_writable(adata):
    obs_index_name = adata.obs.index.name if isinstance(adata.obs.index.name, str) else "obs_names"
    var_index_name = adata.var.index.name if isinstance(adata.var.index.name, str) else "var_names"
    if obs_index_name in adata.obs.columns:
        obs_index_name = "obs_names"
    if var_index_name in adata.var.columns:
        var_index_name = "var_names"
    adata.obs.index = pd.Index(adata.obs.index.astype(str).astype(object), name=obs_index_name)
    adata.var.index = pd.Index(adata.var.index.astype(str).astype(object), name=var_index_name)

    def sanitize(value):
        if isinstance(value, pd.DataFrame):
            if not isinstance(value.index.name, str):
                value.index.name = "index"
            value.index = pd.Index(value.index.astype(str).astype(object), name=value.index.name)
            value.columns.name = str(value.columns.name) if value.columns.name is not None else None
            for col in value.columns:
                if str(value[col].dtype).startswith("string"):
                    value[col] = value[col].astype(str).astype(object)
        elif isinstance(value, dict):
            for sub_value in value.values():
                sanitize(sub_value)

    sanitize(adata.obs)
    sanitize(adata.var)
    sanitize(adata.uns)


def plot_cluster(
    adata,
    barcode_a_whitelist,
    barcode_b_whitelist,
    output,
    config,
    cb_len,
):
    # extract umi counts and gene counts before normalization
    result = adata.obs[['total_counts', 'n_genes_by_counts']].reset_index()
    result = result.rename(columns={result.columns[0]: 'barcode'})

    n_comps = sct_normalize_and_pca(adata)
    print(f"Number of highly variable genes: {adata.var['highly_variable'].sum()}")
    print('\n')

    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=n_comps, show=False)
    ax = plt.gca()
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_title('Variance ratio')
    set_axis_font_sizes(ax)
    plt.savefig(f'{output}/pca.png', bbox_inches='tight', dpi=300)
    plt.close()

    build_snn_graph(adata, n_neighbors=30, n_pcs=20)
    sc.tl.umap(adata, random_state=RANDOM_STATE)
    sc.tl.leiden(
        adata,
        resolution=0.2,
        adjacency=adata.obsp["connectivities"],
        flavor="igraph",
        n_iterations=2,
        directed=False,
        random_state=RANDOM_STATE,
    )
    print(f"Number of clusters: {len(adata.obs['leiden'].unique())}")
    print(f"Cluster sizes:\n{adata.obs['leiden'].value_counts().sort_index()}")

    cluster_ids = sorted(adata.obs['leiden'].astype(int).unique())
    cluster_categories = [str(cluster_id) for cluster_id in cluster_ids]
    adata.obs['leiden'] = pd.Categorical(
        adata.obs['leiden'].astype(str),
        categories=cluster_categories,
        ordered=True,
    )
    cluster_hex_colors = categorical_colors(len(cluster_ids))
    cluster_colors = {
        cluster_id: np.asarray(to_rgba(cluster_hex_colors[i]))
        for i, cluster_id in enumerate(cluster_ids)
    }
    adata.uns['leiden_colors'] = cluster_hex_colors

    sc.pl.umap(adata, color='leiden', legend_loc='on data', title='UMAP - Clusters', frameon=False, show=False)
    ax = plt.gca()
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_title('UMAP')
    set_axis_font_sizes(ax)
    plt.savefig(f'{output}/umap.png', bbox_inches='tight', dpi=300)
    plt.close()

    result = add_spatial_coordinates(
        result,
        "barcode",
        barcode_a_whitelist,
        cb_len,
        barcode_b_whitelist_path=barcode_b_whitelist,
    )
    result = result.merge(adata.obs[['leiden']], left_on='barcode', right_index=True)
    data = pd.DataFrame(
        {
            'x': result['x'],
            'y': result['y'],
            'umi_count': result['total_counts'],
            'gene_count': result['n_genes_by_counts'],
            'leiden': result['leiden'].astype(int),
        }
    )
    data['color'] = data['leiden'].map(lambda cluster_id: to_hex(cluster_colors[int(cluster_id)]))

    spatial_metadata = result.set_index('barcode').loc[adata.obs_names]
    adata.obs['xbc'] = spatial_metadata['xbc'].to_numpy()
    adata.obs['ybc'] = spatial_metadata['ybc'].to_numpy()
    adata.obs['x'] = spatial_metadata['x'].to_numpy(dtype=int)
    adata.obs['y'] = spatial_metadata['y'].to_numpy(dtype=int)
    adata.obs['color'] = adata.obs['leiden'].map(
        lambda cluster_id: to_hex(cluster_colors[int(cluster_id)])
    )
    adata.obsm['spatial'] = adata.obs[['x', 'y']].to_numpy()

    frame_umap = np.zeros(spatial_frame_shape(config), dtype=np.uint8)

    for _, row in data.iterrows():
        x_idx = int(row['x'])
        y_idx = int(row['y'])
        id = int(row['leiden'])

        x_start = x_idx * (config.length_spot + config.interval)
        y_start = y_idx * (config.length_spot + config.interval)
        x_end = x_start + config.length_spot
        y_end = y_start + config.length_spot
        frame_umap[y_start:y_end, x_start:x_end, :] = (cluster_colors[id] * 255).astype(np.uint8)

    img_umap = Image.fromarray(frame_umap, mode = 'RGBA')
    output_size = spatial_image_size(config)
    img_umap = img_umap.resize(output_size, resample=Image.NEAREST)
    box_width = 10
    width, height = img_umap.size
    left = 0
    top = 0
    right = width
    bottom = height
    draw = ImageDraw.Draw(img_umap)
    draw.rectangle(
        [left, top, right, bottom],
        outline=THRESHOLD_COLOR,
        width=box_width,
    )
    img_umap.save(f'{output}/frame_umap.png')
    mask = np.zeros((output_size[1], output_size[0]), dtype=np.uint8)
    mask = Image.fromarray(mask, mode = 'L')
    mask.save(f'{output}/mask.png')

    fig_leg, ax_leg = plt.subplots(figsize=(2, 4))
    ax_leg.axis('off')
    legend_elements = [
        Patch(
            facecolor=cluster_colors[cluster_id],
            label=f'Cluster {cluster_id}'
        )
        for cluster_id in cluster_ids
    ]
    ax_leg.legend(handles=legend_elements, loc='center', frameon=False)
    fig_leg.savefig(f'{output}/umap_legend.png', bbox_inches='tight', dpi=300)
    plt.close(fig_leg)

    make_h5ad_names_writable(adata)
    if 'gene_name' in adata.var:
        adata.var['gene_name'] = adata.var['gene_name'].astype('string')
    ad.settings.allow_write_nullable_strings = True
    h5ad_path = f'{output}/clustered.h5ad'
    adata.write_h5ad(h5ad_path, convert_strings_to_categoricals=False)
    data.to_csv(f'{output}/data.csv', index=False)

    return f'{output}/data.csv'


def find_gene_directory(file_path):
    root = Path(file_path)
    if (root / "GeneFull").is_dir():
        return root / "GeneFull"
    if (root / "Gene").is_dir():
        return root / "Gene"
    raise ValueError("No GeneFull or Gene folder found in the directory")


def run_mrna_qc(
    config,
    file_path,
    barcode_a_whitelist,
    barcode_b_whitelist,
    cb_len,
    umi_min,
    gene_min,
    min_cells,
):
    validate_barcode_length(cb_len)
    method_path = find_gene_directory(file_path) / "raw"
    adata = load_and_filter_counts(
        method_path, True, umi_min, gene_min, min_cells
    )
    csv_path = plot_cluster(
        adata.copy(),
        barcode_a_whitelist,
        barcode_b_whitelist,
        method_path,
        config,
        cb_len,
    )
    plot_spatial_heatmaps(
        csv_path,
        barcode_a_whitelist,
        barcode_b_whitelist,
        method_path,
        cb_len,
        config.x_spots_number,
        config.y_spots_number,
    )


def main():
    args = parse_args()
    config = SpatialPlotConfig(
        args.x_spots_number,
        args.y_spots_number,
        args.length_spot,
        args.interval,
        args.pixel_length,
    )
    run_mrna_qc(
        config,
        args.file_path,
        args.barcodeA_whitelist,
        args.barcodeB_whitelist,
        args.cb_len,
        args.umi_min,
        args.gene_min,
        args.min_cells,
    )


if __name__ == "__main__":
    main()
