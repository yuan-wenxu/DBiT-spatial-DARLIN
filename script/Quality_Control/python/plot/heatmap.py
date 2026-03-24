import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import argparse

def plot_heatmap(file_path, whitelist_path, output):
    data = pd.read_csv(file_path, header = 0)
    if 'x' in data.columns:
        if 'umi_count' in data.columns:
            pass
        else:
            data = data.groupby(['x', 'y']).agg({'UR': 'nunique', 'reads': 'sum'}).reset_index()
            data['umi_count'] = data['UR']
    else:
        data['xbc'] = data['SR'].str[8:16]
        data['ybc'] = data['SR'].str[:8]
        data['UMI'] = data['UR']
        data['reads'] = data['reads']
        whitelist = [line.strip() for line in open(whitelist_path).readlines()]
        data['x'] = data['xbc'].apply(lambda bc: whitelist.index(bc) if bc in whitelist else -1)
        data['y'] = data['ybc'].apply(lambda bc: whitelist.index(bc) if bc in whitelist else -1)
        data = data.groupby(['x', 'y']).agg({'UR': 'nunique', 'reads': 'sum'}).reset_index()
        data['umi_count'] = data['UR']

    # Complete the missing combination
    full_index = pd.MultiIndex.from_product(
        [range(50), range(50)],
        names=['x', 'y']
        )

    data_read = (
        data
        .set_index(['x', 'y'])
        .reindex(full_index, fill_value=0)
        .reset_index()
        )
    try:
        data_read_pivot = data_read.pivot_table(index='y', columns='x', values='reads', aggfunc='sum', fill_value=0)
        plt.figure(figsize=(5,4))
        ax = sns.heatmap(data_read_pivot, cmap='YlOrRd', annot=False, fmt='g', linewidths=0, xticklabels=False, yticklabels=False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_title(f"Reads counts", fontsize = 14)
        plt.tight_layout()
        plt.savefig(os.path.join(output, f"Reads_counts_heatmap.png"), dpi=600)
    except:
        print("No reads data found in the file")

    try:
        data_umi_pivot = data_read.pivot_table(index='y', columns='x', values='umi_count', aggfunc='sum', fill_value=0)
        plt.figure(figsize=(5,4))
        ax = sns.heatmap(data_umi_pivot, cmap='YlOrRd', annot=False, fmt='g', linewidths=0, xticklabels=False, yticklabels=False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_title(f"UMI counts", fontsize = 14)
        plt.tight_layout()
        plt.savefig(os.path.join(output, f"UMI_counts_heatmap.png"), dpi=600)
    except:
        print("No UMI data found in the file")

    try:
        data_umi_pivot = data_read.pivot_table(index='y', columns='x', values='gene_count', aggfunc='sum', fill_value=0)
        plt.figure(figsize=(5,4))
        ax = sns.heatmap(data_umi_pivot, cmap='YlOrRd', annot=False, fmt='g', linewidths=0, xticklabels=False, yticklabels=False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_title(f"Gene counts", fontsize = 14)
        plt.tight_layout()
        plt.savefig(os.path.join(output, f"Gene_counts_heatmap.png"), dpi=600)
    except:
        print("No gene data found in the file")

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Generate heatmap of reads and UMI counts')
    parser.add_argument('-f', '--file_path', type=str, help='Path to the directory containing the final.csv file')
    parser.add_argument('-w', '--whitelist_path', type=str, help='Path to the whitelist file')
    parser.add_argument('-o', '--output', type=str, help='Path to the output directory')
    args = parser.parse_args()

    plot_heatmap(args.file_path, args.whitelist_path, args.output)