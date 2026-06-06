import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ast
import argparse

def filter(file_path, cutoff):

    cell_file = file_path + '/' + 'cell_num_area.csv'
    output_path = file_path + '/' + 'area.png'
    result_path = file_path + '/' + 'filtered_results.csv'

    file = pd.read_csv(cell_file)
    file['area'] = file['area'].apply(ast.literal_eval)
    all_areas = np.concatenate([np.array(sub) for sub in file['area'] if len(sub) > 0])
    print(f'Cell number before filtering: {len(all_areas)}')
    print(f"Spots number before filtering: {(file['num_cells'] != 0).sum()}")

    plt.figure(figsize=(5, 4))
    plt.hist(all_areas, bins=100)
    plt.axvline(cutoff, linestyle='--', linewidth=1, color='r', label=f'cutoff = {cutoff}')
    plt.xlabel('Area')
    plt.ylabel('Count')
    plt.title('Cell area distribution')
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')

    file['count'] = [sum(v >= cutoff for v in sublist) for sublist in file['area']]
    print(f'Cell number after filtering: {file["count"].sum()}')
    print(f"Spots number after filtering: {(file['count'] != 0).sum()}")

    file.to_csv(result_path, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter cells based on area cutoff.')
    parser.add_argument('-f', '--file_path', type=str, help='Path to the directory containing cell_num_area.csv')
    parser.add_argument('-c', '--cutoff', type=float, help='Area cutoff for filtering cells')
    args = parser.parse_args()

    filter(args.file_path, args.cutoff)