import pandas as pd
import argparse
import warnings
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
from matplotlib.patches import Patch
import numpy as np
from PIL import Image, ImageDraw
from scipy import sparse
from sklearn.neighbors import NearestNeighbors

RANDOM_STATE = 42

class PlotConfig:
    def __init__(self, x_spots_number, y_spots_number, length_spot, interval, pixel_length):
        self.x_spots_number = x_spots_number
        self.y_spots_number = y_spots_number
        self.length_spot = length_spot
        self.interval = interval
        self.pixel_length = pixel_length


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


def plot_cluster(adata, whitelist_path, output, config):

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
    ax.set_title('Variance ratio', fontsize=16)
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
    cmap = plt.get_cmap('tab20', len(cluster_ids))
    colors = cmap(range(len(cluster_ids)))
    cluster_colors = {
        cluster_id: colors[i]
        for i, cluster_id in enumerate(cluster_ids)
    }
    adata.uns['leiden_colors'] = [to_hex(colors[i]) for i in range(len(cluster_ids))]

    sc.pl.umap(adata, color='leiden', legend_loc='on data', title='UMAP - Clusters', frameon=False, show=False)
    ax = plt.gca()
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_title('UMAP', fontsize=16)
    plt.savefig(f'{output}/umap.png', bbox_inches='tight', dpi=300)
    plt.close()

    result['xbc'] = result['barcode'].str[8:16]
    result['ybc'] = result['barcode'].str[:8]
    whitelist = [line.strip() for line in open(whitelist_path).readlines()]
    result['x'] = result['xbc'].apply(lambda bc: whitelist.index(bc) if bc in whitelist else -1)
    result['y'] = result['ybc'].apply(lambda bc: whitelist.index(bc) if bc in whitelist else -1)
    result = result.merge(adata.obs[['leiden']], left_on='barcode', right_index=True)
    data = pd.DataFrame(
        {'x': result['x'], 
         'y': result['y'], 
         'umi_count': result['total_counts'], 
         'gene_count': result['n_genes_by_counts'], 
         'leiden': result['leiden'].astype(int)}
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

    frame_umap = np.zeros((int(config.x_spots_number * config.length_spot + (config.x_spots_number-1) * config.interval),
                           int(config.y_spots_number * config.length_spot + (config.y_spots_number-1) * config.interval), 4),
                           dtype=np.uint8)

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
    img_umap = img_umap.resize((int((config.x_spots_number * config.length_spot + (config.x_spots_number - 1) * config.interval)/config.pixel_length),
                                int((config.y_spots_number * config.length_spot + (config.y_spots_number - 1) * config.interval)/config.pixel_length)),
                                resample = Image.NEAREST)
    box_width = 10
    width, height = img_umap.size
    left = 0
    top = 0
    right = width
    bottom = height
    draw = ImageDraw.Draw(img_umap)
    draw.rectangle([left, top, right, bottom], outline="red", width=box_width)
    img_umap.save(f'{output}/frame_umap.png')
    mask = np.zeros((int((config.x_spots_number * config.length_spot + (config.x_spots_number - 1) * config.interval)/config.pixel_length), 
                     int((config.y_spots_number * config.length_spot + (config.y_spots_number - 1) * config.interval)/config.pixel_length)), 
                     dtype=np.uint8)
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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot cluster')
    parser.add_argument('-f', '--file_path', type=str, help='Path to directory containing GeneFull or Gene folder', required=True)
    parser.add_argument('-w', '--whitelist_path', type=str, help='Path to whitelist file', required=True)
    parser.add_argument('-o', '--output_path', type=str, help='Path to output directory', required=True)
    parser.add_argument('--x_spots_number', type=int, default=50, help='Number of spots in x direction')
    parser.add_argument('--y_spots_number', type=int, default=50, help='Number of spots in y direction')
    parser.add_argument('--length_spot', type=int, default=20, help='Length of each spot in pixels')
    parser.add_argument('--interval', type=int, default=20, help='Interval between spots in pixels')
    parser.add_argument('--pixel_length', type=float, default=0.294, help='Length of each pixel in microns')
    args = parser.parse_args()

    config = PlotConfig(args.x_spots_number, args.y_spots_number, args.length_spot, args.interval, args.pixel_length)
    data = plot_cluster(args.file_path, args.whitelist_path, args.output_path, config)
