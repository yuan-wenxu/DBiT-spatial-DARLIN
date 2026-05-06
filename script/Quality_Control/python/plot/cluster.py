import pandas as pd
import argparse
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
from PIL import Image, ImageDraw

class PlotConfig:
    def __init__(self, x_spots_number, y_spots_number, length_spot, interval, pixel_length):
        self.x_spots_number = x_spots_number
        self.y_spots_number = y_spots_number
        self.length_spot = length_spot
        self.interval = interval
        self.pixel_length = pixel_length

def plot_cluster(adata, whitelist_path, output, config):

    # extract umi counts and gene counts before normalization
    result = adata.obs[['total_counts', 'n_genes_by_counts']].reset_index()
    result = result.rename(columns={0: 'barcode'})

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    print(f"Number of highly variable genes: {adata.var['highly_variable'].sum()}")
    print('\n')

    sc.tl.pca(adata, svd_solver='arpack', n_comps=50)
    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50, show=False)
    ax = plt.gca()
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_title('Variance ratio', fontsize=16)
    plt.savefig(f'{output}/pca.png', bbox_inches='tight', dpi=600)
    plt.close()

    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.5)
    print(f"Number of clusters: {len(adata.obs['leiden'].unique())}")
    print(f"Cluster sizes:\n{adata.obs['leiden'].value_counts().sort_index()}")

    sc.pl.umap(adata, color='leiden', legend_loc='on data', title='UMAP - Clusters', frameon=False, show=False)
    ax = plt.gca()
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_title('UMAP', fontsize=16)
    plt.savefig(f'{output}/umap.png', bbox_inches='tight', dpi=600)
    plt.close()

    clusters = len(adata.obs['leiden'].unique())
    cmap = plt.get_cmap('tab20', clusters)
    colors = cmap(range(clusters))
   
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
        frame_umap[y_start:y_end, x_start:x_end, :] = (colors[id] * 255).astype(np.uint8)

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
            facecolor=colors[i],
            label=f'Cluster {i}'
        )
        for i in range(clusters)
    ]
    ax_leg.legend(handles=legend_elements, loc='center', frameon=False)
    fig_leg.savefig(f'{output}/umap_legend.png', bbox_inches='tight', dpi=600)
    plt.close(fig_leg)

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