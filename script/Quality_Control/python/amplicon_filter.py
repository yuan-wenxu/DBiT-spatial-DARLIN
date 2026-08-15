import argparse
from pathlib import Path

import pandas as pd

from utils import (
    ScatterConfig,
    SpatialPlotConfig,
    add_spatial_coordinates,
    load_cell_numbers,
    merge_tissue_spots,
    plot_scatter,
    plot_spatial_frames,
    validate_barcode_length,
)

DARLIN_LOCI = ("CA", "RA", "TA")


def filter_amplicon_by_tissue(
    cell_number_file,
    darlin_path,
    umi_config,
    whitelist_path,
    plot_config,
    cb_len,
):
    validate_barcode_length(cb_len)
    cell_number = load_cell_numbers(cell_number_file)
    darlin_path = Path(darlin_path)
    for locus in DARLIN_LOCI:
        locus_path = darlin_path / locus
        darlin_file = locus_path / "final.csv"
        if darlin_file.exists():
            darlin_data = pd.read_csv(darlin_file)
            darlin_data = add_spatial_coordinates(
                darlin_data, "SR", whitelist_path, cb_len
            )
            merge_data = merge_tissue_spots(darlin_data, cell_number)
            merge_data.to_csv(locus_path / "tissuefiltered.csv", index=False)

            umi_data = merge_data[['x', 'y', 'count', 'UR']]
            umi_data = umi_data.groupby(['x', 'y']).agg({'count': 'first', 'UR': 'nunique'}).reset_index()
            umi_data['umi_count'] = umi_data['UR']
            umi_with_cells = umi_data[umi_data['count'] > 0]
            if not umi_with_cells.empty:
                plot_scatter(umi_with_cells['count'], umi_with_cells['umi_count'], locus_path, umi_config)
            plot_spatial_frames(umi_data, locus_path, plot_config)

            darlin_data = merge_data[['x', 'y', 'count', 'n_LR']]
            darlin_with_cells = darlin_data[darlin_data['count'] > 0]
            if not darlin_with_cells.empty:
                plot_scatter(darlin_with_cells['count'], darlin_with_cells['n_LR'], locus_path, ScatterConfig('Number of cells', 'Number of lineage barcodes', 'Lineage_barcode_distribution', False, True, False, False))

            print(locus)
            print(f'Spots number: {len(umi_data)}')
            print(f'UMI number: {umi_data["umi_count"].sum()}')
            if 'LR' in merge_data.columns:
                print(f'Lineage barcode number: {merge_data["LR"].nunique()}')
            print('\n')
        else:
            print(f'{darlin_file} does not exist.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot amplicon data after filtering spots outside tissue.')
    parser.add_argument('-c', '--cell_number_file', type=str, help='cell number file')
    parser.add_argument('-d', '--darlin_path', type=str, help='data path')
    parser.add_argument('-w', '--whitelist_path', type=str, help='whitelist file')
    parser.add_argument('--cb-len', type=int, required=True, help='Total concatenated cell-barcode length in bp')
    parser.add_argument('--x_spots_number', type=int, default=50, help='Number of spots in x direction')
    parser.add_argument('--y_spots_number', type=int, default=50, help='Number of spots in y direction')
    parser.add_argument('--length_spot', type=int, default=20, help='Length of each spot in pixels')
    parser.add_argument('--interval', type=int, default=20, help='Interval between spots in pixels')
    parser.add_argument('--pixel_length', type=float, default=0.294, help='Length of each pixel in microns')
    args = parser.parse_args()

    cell_number_file = args.cell_number_file
    darlin_path = args.darlin_path
    whitelist_path = args.whitelist_path
    x_spots_number = args.x_spots_number
    y_spots_number = args.y_spots_number
    length_spot = args.length_spot
    interval = args.interval
    pixel_length = args.pixel_length

    frame_config = SpatialPlotConfig(x_spots_number, y_spots_number, length_spot, interval, pixel_length)

    umi_config = ScatterConfig('Number of cells', 'Number of UMIs', 'UMI_distribution', False, True, False, False)
    filter_amplicon_by_tissue(
        cell_number_file,
        darlin_path,
        umi_config,
        whitelist_path,
        frame_config,
        args.cb_len,
    )
