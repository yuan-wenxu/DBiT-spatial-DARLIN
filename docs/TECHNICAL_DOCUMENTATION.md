# DBiT-spatial-DARLIN Technical Documentation

This document provides a detailed technical description of the DBiT-spatial-DARLIN pipeline.

- Step scripts:
  - `dbit_mrna.sh` 
  - `image.sh`
  - `dbit_amplicon.sh`
  - `plot_cell_filtered.sh`

---

## 0. Pipeline Overview

This project is used for quality control of sequencing data from DBiT.

There are four shell scripts to process data from three modalities.

`plot_cell_filtered.sh` is used to integrate data from three modalities.

---

## 1. Transcriptome

1. Extract spatial barcode and UMI based on sequence information.
  
    This is the structure of barcode and UMI.

    ![Barcode and UMI structure](image/barcode.png)


2. Using STARsolo to align transcriptome to genome.
3. Cluster spots and plot the cluster labels according to spatial location.

### 1.1 Transcriptome Clustering

The mRNA QC clustering step is implemented in `script/Quality_Control/python/plot/cluster.py`.
It starts from the STARsolo count matrix after QC filtering and uses the following workflow:

1. Save raw QC metrics before normalization.

   The raw total UMI count (`total_counts`) and detected gene count (`n_genes_by_counts`) are extracted before normalization. These values are written to `data.csv` as `umi_count` and `gene_count`.

2. Normalize and reduce dimensionality with an SCT-like Pearson residual workflow.

   The pipeline uses `scanpy.experimental.pp.recipe_pearson_residuals` with:

   - `theta = 100`
   - `n_top_genes = 3000`, capped by the number of available genes
   - `n_comps = 50`, capped by the number of available spots and genes
   - `random_state = 42`

   This is a Pearson residual normalization workflow similar in purpose to SCT normalization, but it is implemented through Scanpy rather than directly through Seurat SCT.

3. Build an SNN graph from PCA coordinates.

   The SNN graph is constructed from `adata.obsm["X_pca"]` using:

   - `n_neighbors = 10`
   - `n_pcs = 20`, capped by the number of available PCs
   - Euclidean nearest-neighbor search
   - Jaccard-style SNN edge weights, calculated as shared neighbors divided by total neighbors

   The resulting graph is stored in:

   - `adata.obsp["connectivities"]`
   - `adata.obsp["distances"]`
   - `adata.uns["neighbors"]`

4. Run UMAP and Leiden clustering.

   UMAP is run on the SNN graph with `random_state = 42`.

   Leiden clustering is run on `adata.obsp["connectivities"]` with:

   - `resolution = 0.2`
   - `flavor = "igraph"`
   - `n_iterations = 2`
   - `directed = False`
   - `random_state = 42`

5. Save clustering outputs.

   The clustering step writes:

   - `pca.png`: PCA variance ratio plot
   - `umap.png`: UMAP colored by Leiden cluster
   - `frame_umap.png`: spatial cluster map
   - `umap_legend.png`: cluster color legend
   - `data.csv`: spot-level QC and cluster table, including `x`, `y`, `umi_count`, `gene_count`, `leiden`, and `color`
   - `clustered.h5ad`: AnnData object with normalized representation, PCA, SNN graph, UMAP, Leiden cluster labels, and cluster colors

## 2. Image

1. Extract the sampling area based on the transcriptome plot results.
2. Cut the image according to the actual sampling area and use stardist to predict the number of cells.

## 3. Amplicon

1. Extract the DARLIN sequence. (optional)
2. Extract spatial barcode and UMI based on sequence information.
3. Correct DARLIN sequence.
