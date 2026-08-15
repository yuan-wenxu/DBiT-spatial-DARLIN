from pathlib import Path

import anndata as ad
import pandas as pd
import numpy as np
import argparse

from utils import (
    ScatterConfig,
    SpatialPlotConfig,
    load_cell_numbers,
    merge_tissue_spots,
    plot_scatter,
    plot_spatial_frames,
)

METHODS = ("raw",)


def filter_clustered_h5ad(h5ad_path, cell_number):
    h5ad_path = Path(h5ad_path)
    if not h5ad_path.exists():
        print(f'Warning: clustered.h5ad not found at {h5ad_path}; skipping tissue filtering.')
        return

    adata = ad.read_h5ad(h5ad_path)
    required = {'x', 'y'}
    missing = required.difference(adata.obs.columns)
    if missing:
        raise ValueError(
            f'clustered.h5ad obs is missing coordinate columns: {sorted(missing)}'
        )

    tissue_coordinates = pd.MultiIndex.from_frame(
        cell_number.loc[cell_number['in_tissue'], ['x', 'y']]
    )
    h5ad_coordinates = pd.MultiIndex.from_frame(adata.obs[['x', 'y']])
    keep = h5ad_coordinates.isin(tissue_coordinates)
    filtered = adata[keep].copy()

    output_path = h5ad_path.with_name(f'{h5ad_path.stem}.tissuefiltered.h5ad')
    filtered.write_h5ad(output_path)

    print(
        f'Tissue-filtered clustered.h5ad: {adata.n_obs} -> {filtered.n_obs} '
        f'spots; saved to {output_path}'
    )


def plot_filtered(cell_number_file, umi_gene, umi_config, gene_config, frame_config):
    cell_number = load_cell_numbers(cell_number_file)
    for method in METHODS:
        method_path = Path(umi_gene) / method
        data_path = method_path / "data.csv"
        data = pd.read_csv(data_path)
        merge_data = merge_tissue_spots(data, cell_number)
        merge_data.to_csv(method_path / "data_tissuefiltered.csv", index=False)
        filter_clustered_h5ad(method_path / "clustered.h5ad", cell_number)
        print(data_path)
        print(f'Tissue-filtered spots: {len(merge_data)}')
        print(f'Total UMI: {np.sum(merge_data["umi_count"])}')
        print(f'Total Gene: {np.sum(merge_data["gene_count"])}')
        print(f"Mean UMI: {np.mean(merge_data['umi_count'])}")
        print(f"Median UMI: {np.median(merge_data['umi_count'])}")
        print(f"Mean Gene: {np.mean(merge_data['gene_count'])}")
        print(f"Median Gene: {np.median(merge_data['gene_count'])}")
        print(f'Number of cells: {np.sum(merge_data["count"])}')
        spots_with_cells = merge_data[merge_data['count'] > 0]
        if not spots_with_cells.empty:
            print(f"Mean UMI per cell: {np.mean(spots_with_cells['umi_count']/spots_with_cells['count'])}")
            print(f"Median UMI per cell: {np.median(spots_with_cells['umi_count']/spots_with_cells['count'])}")
            plot_scatter(spots_with_cells['count'], spots_with_cells['umi_count'], method_path, umi_config)
            plot_scatter(spots_with_cells['count'], spots_with_cells['gene_count'], method_path, gene_config)
        print('\n')

        plot_spatial_frames(merge_data, method_path, frame_config)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot filtered results')
    parser.add_argument('-c', '--cell_number_file', type=str, help='cell number file')
    parser.add_argument('-d', '--data_path', type=str, help='data path')
    parser.add_argument('--x_spots_number', type=int, default=50, help='Number of spots in x direction')
    parser.add_argument('--y_spots_number', type=int, default=50, help='Number of spots in y direction')
    parser.add_argument('--length_spot', type=int, default=20, help='Length of each spot in pixels')
    parser.add_argument('--interval', type=int, default=20, help='Interval between spots in pixels')
    parser.add_argument('--pixel_length', type=float, default=0.294, help='Length of each pixel in microns')
    args = parser.parse_args()

    cell_number_file = args.cell_number_file
    data_path = args.data_path
    x_spots_number = args.x_spots_number
    y_spots_number = args.y_spots_number
    length_spot = args.length_spot
    interval = args.interval
    pixel_length = args.pixel_length

    frame_config = SpatialPlotConfig(x_spots_number, y_spots_number, length_spot, interval, pixel_length)

    umi_config = ScatterConfig('Number of cells', 'Number of UMIs', 'UMI_distribution', False, True, False, False)
    gene_config = ScatterConfig('Number of cells', 'Number of genes', 'Gene_distribution', False, True, False, False)

    plot_filtered(cell_number_file, data_path, umi_config, gene_config, frame_config)
