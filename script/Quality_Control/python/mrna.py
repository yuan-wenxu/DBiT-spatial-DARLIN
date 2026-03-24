from plot import PlotConfig
from plot import plot_cluster
from plot import plot_violin
from plot import plot_heatmap
from utils import extract
import os
import argparse

darlin = ['CA', 'RA', 'TA']
method = ['raw', 'filtered']

def mrna(config, file_path, whitelist_path, umi_min, gene_min, min_cells):
    if os.path.isdir(file_path + '/' + 'GeneFull' ):
        data_path = file_path + '/' + 'GeneFull' 
    elif os.path.isdir(file_path + '/' + 'Gene' ):
        data_path = file_path + '/' + 'Gene' 
    else:
        raise ValueError('No GeneFull or Gene folder found in the directory')
    for m in method:
        if m == 'raw':
            filter_mode = True
        else:
            filter_mode = False
        m_path = data_path + '/' + m
        adata = plot_violin(m_path, filter_mode, umi_min, gene_min, min_cells)
        bdata = adata.copy()
        csv_path = plot_cluster(bdata, whitelist_path, m_path, config)
        plot_heatmap(csv_path, whitelist_path, m_path)
        for d in darlin:
            try:
                darlin_file = extract(adata, m_path, whitelist_path, d)
                plot_heatmap(darlin_file, whitelist_path, f'{m_path}/{d}')
            except:
                print(f'No {d} data found')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot mRNA data')
    parser.add_argument('-f', '--file_path', type=str, help='Path to directory containing GeneFull or Gene folder', required=True)
    parser.add_argument('-w', '--whitelist_path', type=str, help='Path to whitelist file', required=True)
    parser.add_argument('-umi_min', '--umi_min', type=int, default=900, help='Minimum UMI count per spot')
    parser.add_argument('-gene_min', '--gene_min', type=int, default=300, help='Minimum gene count per spot')
    parser.add_argument('-min_cells', '--min_cells', type=int, default=3, help='Minimum number of cells per gene')
    parser.add_argument('--x_spots_number', type=int, default=50, help='Number of spots in x direction')
    parser.add_argument('--y_spots_number', type=int, default=50, help='Number of spots in y direction')
    parser.add_argument('--length_spot', type=int, default=20, help='Length of each spot in pixels')
    parser.add_argument('--interval', type=int, default=20, help='Interval between spots in pixels')
    parser.add_argument('--pixel_length', type=float, default=0.294, help='Length of each pixel in microns')
    args = parser.parse_args()
    config = PlotConfig(args.x_spots_number, args.y_spots_number, args.length_spot, args.interval, args.pixel_length)
    mrna(config, args.file_path, args.whitelist_path, args.umi_min, args.gene_min, args.min_cells)