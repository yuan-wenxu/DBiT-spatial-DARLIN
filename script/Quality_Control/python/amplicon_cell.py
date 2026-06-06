import pandas as pd
from plot import plot_scatter, ScatterConfig, plot_frame_filtered, PlotConfig
import argparse
import os

darlin = ['CA', 'RA', 'TA']

def main(cell_number_file, darlin_path, umi_config, whitelist_path, plot_config):
    cell_number = pd.read_csv(cell_number_file, header = 0)
    whitelist = [line.strip() for line in open(whitelist_path).readlines()]
    for d in darlin:
        darlin_file = darlin_path + '/' + d + '/' + 'final.csv'
        if os.path.exists(darlin_file):
            darlin_data = pd.read_csv(darlin_file, header = 0)
            darlin_data['xbc'] = darlin_data['SR'].str[8:16]
            darlin_data['ybc'] = darlin_data['SR'].str[:8]
            darlin_data['x'] = darlin_data['xbc'].apply(lambda bc: whitelist.index(bc) if bc in whitelist else -1)
            darlin_data['y'] = darlin_data['ybc'].apply(lambda bc: whitelist.index(bc) if bc in whitelist else -1)        
            merge_data = (darlin_data.merge(cell_number[['x', 'y', 'count']], on = ['x', 'y']).fillna(0))
            merge_data = merge_data[merge_data['count'] > 0]
            merge_data.to_csv(darlin_path + '/' + d + '/' + 'cellfiltered.csv', index = False)

            umi_data = merge_data[['x', 'y', 'count', 'UR']]
            umi_data = umi_data.groupby(['x', 'y']).agg({'count': 'first', 'UR': 'count'}).reset_index()
            umi_data['umi_count'] = umi_data['UR']
            plot_scatter(umi_data['count'], umi_data['umi_count'], darlin_path + '/' + d, umi_config)
            plot_frame_filtered(umi_data, darlin_path + '/' + d, plot_config)

            print(d)
            print(f'Spots number: {len(umi_data)}')
            print(f'UMI number: {umi_data["umi_count"].sum()}')
            if 'LR' in merge_data.columns:
                print(f'Lineage barcode number: {merge_data["LR"].nunique()}')
            print('\n')
        else:
            print(f'{darlin_file} does not exist.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot amplicon data after filtering out spots without cells.')
    parser.add_argument('-c', '--cell_number_file', type=str, help='cell number file')
    parser.add_argument('-d', '--darlin_path', type=str, help='data path')
    parser.add_argument('-w', '--whitelist_path', type=str, help='whitelist file')
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

    frame_config = PlotConfig(x_spots_number, y_spots_number, length_spot, interval, pixel_length)

    umi_config = ScatterConfig('Number of cells', 'Number of UMIs', 'UMI_distribution', False, True, False, False)
    main(cell_number_file, darlin_path, umi_config, whitelist_path, frame_config)
