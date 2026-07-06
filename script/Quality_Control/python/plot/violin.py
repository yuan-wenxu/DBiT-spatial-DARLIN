from scipy.io import mmread
import os
import pandas as pd
import argparse
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def plot_violin(file_path, filter, umi_min=900, gene_min=300, min_cells=3):
    mtx = mmread(file_path+ '/' + 'matrix.mtx')
    barcodes = pd.read_csv(os.path.join(file_path, "barcodes.tsv"), header=None, sep='\t', dtype=str)
    features = pd.read_csv(os.path.join(file_path, "features.tsv"), header=None, sep='\t', dtype=str)
    adata = sc.AnnData(X = mtx.T.tocsr())
    adata.obs_names = barcodes.iloc[:, 0].astype(str).to_numpy()
    adata.var_names = features.iloc[:, 0].astype(str).to_numpy()
    adata.obs_names.name = "barcode"
    adata.var_names.name = "gene"
    sc.pp.calculate_qc_metrics(adata, inplace=True)

    # Violin plot
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    sns.violinplot(y=adata.obs['n_genes_by_counts'], ax=axes[0], color = '#b2df8a',inner="box", width = 0.8)
    axes[0].set_title('Genes per Spot', fontsize = 12)
    axes[0].set_ylabel('Number of Genes', fontsize = 10)
    axes[0].tick_params(axis='y', labelsize=10)
    sns.violinplot(y=adata.obs['total_counts'], ax=axes[1], color='#ffffb3', inner="box", width = 0.8)
    axes[1].set_title('UMI per Spot', fontsize = 12)
    axes[1].set_ylabel('Number of UMIs', fontsize = 10)
    axes[1].tick_params(axis='y', labelsize=10)
    plt.subplots_adjust(wspace=0.4)
    plt.tight_layout()
    plt.savefig(f'{file_path}/violin.png', bbox_inches="tight", dpi=300)
    plt.close()

    print(f'Before filtering: {adata.n_obs} spots, {adata.n_vars} genes')

    if filter:
        # Filter cells based on QC metrics
        gene_counts = adata.obs['n_genes_by_counts'].to_numpy()
        gene_count_p95 = np.percentile(gene_counts, 95)
        plt.figure(figsize=(5, 4))
        plt.hist(gene_counts[gene_counts <= gene_count_p95], bins=100, color='#b2df8a')
        plt.axvline(gene_min, color='r', linestyle='--', label=f'Gene Minimum: {gene_min}')
        plt.xlim(0, gene_count_p95)
        plt.legend(loc='upper right')
        plt.title(
            f'Genes per Spot\nUpper 5% omitted (P95 = {gene_count_p95:,.0f})',
            fontsize=12,
        )
        plt.xlabel('Number of Genes', fontsize=10)
        plt.ylabel('Frequency', fontsize=10)
        plt.tick_params(axis='both', labelsize=10)
        plt.tight_layout()
        plt.savefig(f'{file_path}/gene_counts_hist.png', bbox_inches="tight", dpi=300)
        plt.close()

        umi_counts = adata.obs['total_counts'].to_numpy()
        umi_count_p95 = np.percentile(umi_counts, 95)
        plt.figure(figsize=(5, 4))
        plt.hist(umi_counts[umi_counts <= umi_count_p95], bins=100, color='#b2df8a')
        plt.axvline(umi_min, color='r', linestyle='--', label=f'UMI Minimum: {umi_min}')
        plt.xlim(0, umi_count_p95)
        plt.legend(loc='upper right')
        plt.title(
            f'UMI per Spot\nUpper 5% omitted (P95 = {umi_count_p95:,.0f})',
            fontsize=12,
        )
        plt.xlabel('Number of UMIs', fontsize=10)
        plt.ylabel('Frequency', fontsize=10)
        plt.tick_params(axis='both', labelsize=10)
        plt.tight_layout()
        plt.savefig(f'{file_path}/umi_counts_hist.png', bbox_inches="tight", dpi=300)
        plt.close()

        spots_per_gene = np.array((adata.X > 0).sum(axis=0)).flatten()
        full_bin_width = max(1, int(np.ceil(adata.n_obs / 100)))
        full_bins = np.arange(0, adata.n_obs + full_bin_width, full_bin_width)
        plt.figure(figsize=(5, 4))
        plt.hist(spots_per_gene, bins=full_bins, color='#b2df8a')
        plt.axvline(min_cells, color='r', linestyle='--', label=f'Minimum Cells: {min_cells}')
        plt.legend(loc='upper right')
        plt.title('Spots per Gene', fontsize=12)
        plt.xlabel('Number of Spots', fontsize=10)
        plt.ylabel('Frequency', fontsize=10)
        plt.tick_params(axis='both', labelsize=10)
        plt.tight_layout()
        plt.savefig(f'{file_path}/spots_per_gene_hist.png', bbox_inches="tight", dpi=300)
        plt.close()

        spots_per_gene = np.array((adata.X > 0).sum(axis=0)).flatten()
        zoom_bins = np.arange(0.5, 21.5, 1)
        plt.figure(figsize=(5, 4))
        plt.hist(spots_per_gene, bins=zoom_bins, color='#b2df8a')
        plt.axvline(min_cells, color='r', linestyle='--', label=f'Minimum Cells: {min_cells}')
        plt.title('Spots per Gene', fontsize=12)
        plt.xlabel('Number of Spots', fontsize=10)
        plt.ylabel('Frequency', fontsize=10)
        plt.tick_params(axis='both', labelsize=10)
        plt.xlim(0.5, 20.5)
        plt.xticks(range(1, 21))
        plt.legend(loc='upper right')
        plt.tight_layout()
        plt.savefig(f'{file_path}/spots_per_gene_hist_small.png', bbox_inches="tight", dpi=300)
        plt.close()

        adata = adata[(adata.obs['total_counts'] >= umi_min) &
              (adata.obs['n_genes_by_counts'] >= gene_min), :].copy()
        sc.pp.filter_genes(adata, min_cells=min_cells)

        print(f'After filtering: {adata.n_obs} spots, {adata.n_vars} genes')

        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        sns.violinplot(y=adata.obs['n_genes_by_counts'], ax=axes[0], color = '#b2df8a',inner="box", width = 0.8)
        axes[0].set_title('Genes per Spot', fontsize = 12)
        axes[0].set_ylabel('Number of Genes', fontsize = 10)
        axes[0].tick_params(axis='y', labelsize=10)
        sns.violinplot(y=adata.obs['total_counts'], ax=axes[1], color='#ffffb3', inner="box", width = 0.8)
        axes[1].set_title('UMI per Spot', fontsize = 12)
        axes[1].set_ylabel('Number of UMIs', fontsize = 10)
        axes[1].tick_params(axis='y', labelsize=10)
        plt.subplots_adjust(wspace=0.4)
        plt.tight_layout()
        plt.savefig(f'{file_path}/violin_filtered_.png', bbox_inches="tight", dpi=300)
        plt.close()
        
    return adata

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Draw violin plot based on matrix.mtx')
    parser.add_argument('-f', '--file_path', type=str, help='Path to the directory of matrix.mtx')
    parser.add_argument('-filter', '--filter', type=bool, help='Filter cells based on QC metrics')
    parser.add_argument('-umi_min', '--umi_min', type=int, default=900, help='Minimum UMI count per spot')
    parser.add_argument('-gene_min', '--gene_min', type=int, default=300, help='Minimum gene count per spot')
    parser.add_argument('-min_cells', '--min_cells', type=int, default=3, help='Minimum number of cells per gene')
    args = parser.parse_args()

    adata = plot_violin(args.file_path, args.filter, args.umi_min, args.gene_min, args.min_cells)
