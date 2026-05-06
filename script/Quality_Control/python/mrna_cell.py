import pandas as pd
import numpy as np
from plot import plot_scatter, ScatterConfig, plot_frame_filtered, PlotConfig
import argparse

method = ['raw', 'filtered']
darlin = ['CA', 'RA', 'TA']

def plot_filtered(cell_number_file, umi_gene, umi_config, gene_config, frame_config):
    cell_number = pd.read_csv(cell_number_file, header = 0)
    #cell_number['y'] = 49 - cell_number['y']
    for m in method:
        m_path = umi_gene + '/' + m
        data_path = m_path + '/' + 'data.csv'
        data = pd.read_csv(data_path, header = 0)
        merge_data = (data.merge(cell_number[['x', 'y', 'count']], on = ['x', 'y']).fillna(0))
        merge_data = merge_data[merge_data['count'] > 0]
        merge_data.to_csv(m_path + '/' + 'data_cellfiltered.csv', index = False)
        print(data_path)
        print(f'Total number of spots: {len(merge_data)}')
        print(f'Total UMI: {np.sum(merge_data["umi_count"])}')
        print(f'Total Gene: {np.sum(merge_data["gene_count"])}')
        print(f"Mean UMI: {np.mean(merge_data['umi_count'])}")
        print(f"Median UMI: {np.median(merge_data['umi_count'])}")
        print(f"Mean Gene: {np.mean(merge_data['gene_count'])}")
        print(f"Median Gene: {np.median(merge_data['gene_count'])}")
        print(f'Number of cells: {np.sum(merge_data["count"])}')
        print(f"Mean UMI per cell: {np.mean(merge_data['umi_count']/merge_data['count'])}")
        print(f"Median UMI per cell: {np.median(merge_data['umi_count']/merge_data['count'])}")
        print(f"Mean Gene per cell: {np.mean(merge_data['gene_count']/merge_data['count'])}")
        print(f"Median Gene per cell: {np.median(merge_data['gene_count']/merge_data['count'])}")
        print('\n')
        plot_scatter(merge_data['count'], merge_data['umi_count'], m_path, umi_config)
        plot_scatter(merge_data['count'], merge_data['gene_count'], m_path, gene_config)

        plot_frame_filtered(merge_data, m_path, frame_config)
        
        for d in darlin:
            d_path = m_path + '/' + d + '/' + f'{d}.csv'
            darlin_data = pd.read_csv(d_path, header = 0)
            darlin_merge = (darlin_data.merge(cell_number[['x', 'y', 'count']], on = ['x', 'y']).fillna(0))
            darlin_merge = darlin_merge[darlin_merge['count'] > 0]
            print(f"Spots number for {d}: {len(darlin_merge)}")
            print(f"Total UMI for {d}: {np.sum(darlin_merge['umi_count'])}")
            print(f'Darlin UMI per spot for {d}: {np.sum(darlin_merge["umi_count"])/len(merge_data)}')
            plot_scatter(darlin_merge['count'], darlin_merge['umi_count'], m_path + '/' + d, umi_config)
            darlin_merge.to_csv(m_path + '/' + d + '/' + f'{d}_cellfiltered.csv', index = False)
        print('\n')

        

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

    frame_config = PlotConfig(x_spots_number, y_spots_number, length_spot, interval, pixel_length)

    umi_config = ScatterConfig('Number of cells', 'Number of UMIs', 'UMI_distribution', False, True, False, False)
    gene_config = ScatterConfig('Number of cells', 'Number of genes', 'Gene_distribution', False, True, False, False)

    plot_filtered(cell_number_file, data_path, umi_config, gene_config, frame_config)