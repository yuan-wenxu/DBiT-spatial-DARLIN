# DBiT-spatial-DARLIN Technical Documentation

This document describes the implementation details of the DBiT-spatial-DARLIN processing workflow. For runnable examples and a shorter overview, see the repository [README](../README.md).

## 1. Pipeline Scope

The repository contains two related workflows:

1. Quality control for DBiT spatial data:
   - transcriptome FASTQ processing and spatial mRNA QC
   - registered image splitting and StarDist cell counting
   - DARLIN amplicon FASTQ processing and spatial clone-call QC
   - tissue-filtered visualization and optional image merging
2. Clone analysis:
   - filter LR sequences against allele-bank definitions
   - plot the top LR clones on the mRNA Leiden cluster background

Primary shell entry points:

```text
script/dbit.sh
script/Quality_Control/mrna.sh
script/Quality_Control/image.sh
script/Quality_Control/amplicon.sh
script/Quality_Control/plot.sh
script/Clone_Analysis/top_lr_pipeline.sh
```

`dbit.sh` is the single-step local/SLURM launcher. For `mrna`, `amplicon`, and
`image`, it receives an input path through `--input` and appends that input plus
derived result paths to the per-dataset config. Later runs may omit `--input`
and reuse the path stored for that step. The final `plot` step runs after image,
reads those accumulated paths, and therefore needs no separate input. The
`clone` step runs `top_lr_pipeline.sh`, defaults to the stored `amp_dir` and
`mrna_dir`, and stores command-line overrides for reuse. `--chip`
accepts `50-50`, `50-20`, or
`100-20`. Chip dimensions and whitelist selection are defined only in that
launcher and exported to the selected worker job. All shell entry points run
Python modules through `pixi run -e <env>`.

The selected chip name is also appended to the per-dataset config. Later steps
may omit `--chip` and reuse that value; a new command-line value is appended as
the latest selection. Grid dimensions remain centralized in `dbit.sh`.

The mRNA step optionally accepts `--genome-dir`, `--umi-min`, `--gene-min`, and
`--min-cell`; provided values are appended to the per-dataset config before
submission and are rejected for other steps. The genome directory override is
stored as an absolute path. The amplicon step similarly accepts
`--initial-reads-cutoff`, `--major-fraction-threshold-molecule`,
`--reads-cutoff`, and `--slope-cutoff`, all restricted to that step. Image runs
require `--orientation` and `--swap-xy` only when the corresponding values are
not already stored; command-line values are appended so later image and plot
runs use the same transformation. These options are rejected for all non-image
steps, and plot fails if image has not stored them.

Clone analysis is also available through `dbit clone`. Like the QC plot step, it
sources the config file and reads all parameters from it, including `clone_labels`
and `clone_top_n`. The launcher uses the `default`
Pixi environment and the clone SLURM settings from the config.

## 2. Shared Concepts

### Spatial Coordinates

Spatial spot coordinates are stored as integer `x` and `y` columns. Most spatial plotting code uses image-like coordinates:

- `x` increases left to right.
- `y` increases top to bottom.
- Matplotlib plots call `invert_yaxis()` when needed to preserve this image-coordinate interpretation.

### Barcode Structure

Transcriptome and amplicon preprocessing both extract a 16 bp spatial barcode and a 10 bp UMI. The 16 bp spatial barcode is built from two 8 bp components.

![Barcode and UMI structure](image/barcode.png)

### Orientation Parameters

`image.sh` and `plot.sh` share orientation controls:

The shared QC config sets `orientation` to `normal`, `horizontal`, `vertical`,
or `rotate`; `swap_xy=True` additionally swaps the coordinate axes.

These parameters are documented in detail in [ORIENTATION.md](ORIENTATION.md). The clone-analysis pipeline reads `orientation` and `swap_xy` from the same shared config, matching the QC pipeline conventions.

## 3. Transcriptome Workflow

Entry point:

```text
script/Quality_Control/mrna.sh
```

### 3.1 Preprocessing

The shell script locates transcriptome FASTQ pairs from the input path passed
by `dbit.sh`, then calls the preprocessing Python code under:

```text
script/Quality_Control/python/preprocessing/
```

Main operations:

1. Match linker sequences.
2. Extract barcode and UMI sequence.
3. Optionally correct barcode components against the whitelist.
4. Write barcode/UMI FASTQ files for downstream STARsolo processing.

Important parameters:

- `linker1`, `linker2`: expected linker sequences.
- `mm_rate`: mismatch rate for linker matching.
- `bc_max_dist`: maximum barcode correction distance.
- `gzip_output`: whether preprocessing output is gzipped immediately.
- `gzip_after_preprocess`: whether uncompressed preprocessing output is compressed after extraction.

### 3.2 STARsolo Alignment

`mrna.sh` runs STAR with `GeneFull` solo features. Relevant STARsolo settings include:

```text
--soloType CB_UMI_Simple
--soloCBstart 1
--soloCBlen 16
--soloUMIstart 17
--soloUMIlen 10
--soloCBwhitelist None
--soloCellFilter EmptyDrops_CR
--soloFeatures GeneFull
```

The default barcode and UMI positions match the preprocessing FASTQ layout:

- cell barcode: bases 1-16
- UMI: bases 17-26

STAR outputs are written under:

```text
<output_path>/results/<sample_name>/Solo.out/GeneFull/
```

### 3.3 mRNA QC and Clustering

After STARsolo, `mrna.sh` calls:

```text
script/Quality_Control/python/mrna.py
```

This runs QC, filtering, plotting, and clustering. The clustering implementation is in:

```text
script/Quality_Control/python/plot/cluster.py
```

The clustering workflow:

1. Load the STARsolo count matrix.
2. Save raw QC metrics before normalization.
3. Filter spots using UMI and gene thresholds.
4. Normalize with a Scanpy Pearson-residual workflow.
5. Run PCA.
6. Build an SNN graph from PCA coordinates.
7. Run UMAP and Leiden clustering.
8. Write spatial cluster plots and tabular outputs.

Key defaults:

- `--umi_min`: `900`
- `--gene_min`: `300`
- `--min_cells`: `3`
- Pearson residual `theta`: `100`
- Pearson residual `n_top_genes`: `3000`, capped by available genes
- PCA `n_comps`: `50`, capped by available spots and genes
- SNN `n_neighbors`: `10`
- SNN `n_pcs`: `20`, capped by available PCs
- Leiden `resolution`: `0.2`
- Leiden `random_state`: `42`

Important mRNA outputs:

```text
Solo.out/GeneFull/raw/
├── data.csv
├── clustered.h5ad
├── pca.png
├── umap.png
├── frame_umap.png
├── umap_legend.png
├── umi_filtered.png
├── gene_filtered.png
└── gene_per_cell_filtered.png
```

`data.csv` contains spot-level information including:

- `x`, `y`
- raw QC metrics such as `umi_count` and `gene_count`
- `leiden`
- `color`

`data_tissuefiltered.csv` is produced later by `plot.sh` after applying the image-derived tissue mask to mRNA spot data.

## 4. Image Workflow

Entry point:

```text
script/Quality_Control/image.sh
```

The image workflow assumes a registered and cropped image, usually named `align.png`.

Main steps:

1. Generate a coarse whole-image tissue mask by treating near-black pixels as
   background, requiring a small amount of local signal density, and removing
   isolated regions outside the main tissue.
2. Pass the mask into the image-splitting loop so each logical DBiT spot is
   numbered, cropped, and classified in one traversal using the same
   orientation mapping.
3. Run StarDist on each generated tile.
4. Write cell count, cell area, and `in_tissue` status for every spot.
5. Filter predicted cells by the configured area cutoff while retaining both
   the per-spot cell-count interface and `in_tissue` column. Downstream plots
   filter spots using `in_tissue` rather than requiring `count > 0`.

Implementation files:

```text
script/Quality_Control/python/stardist_segment.py
script/Quality_Control/python/image_process/tissue_mask.py
script/Quality_Control/python/image_process/split.py
script/Quality_Control/python/image_process/stardist_predict.py
script/Quality_Control/python/cell_filter.py
```

Important outputs:

```text
image/
├── result.png
├── cell_num_area.csv
├── filtered_results.csv
├── tissue_mask.png
├── mask/
├── label/
└── split/
```

`filtered_results.csv` is the image-derived table used by `plot.sh`. It retains
spot coordinates and predicted cell counts and adds the Boolean `in_tissue`
column used for spot filtering. `tissue_mask.png` is a full-resolution binary
mask in the registered image coordinate system.

Tissue-mask thresholds and cleanup parameters are internal defaults; no
additional user configuration is required.

The image workflow uses the `image` Pixi environment because it depends on TensorFlow, StarDist, OpenCV, and related image packages.

## 5. Amplicon Workflow

Entry point:

```text
script/Quality_Control/amplicon.sh
```

The amplicon workflow processes CA, RA, and TA DARLIN amplicon FASTQ files. The
script infers the locus from sample names containing `CA`, `RA`, or `TA`, and
accepts either `sample-CA` or `sample_CA` naming.

### 5.1 Preprocessing

Main operations:

1. Optionally run cutadapt.
2. Extract spatial barcode and UMI.
3. Treat the DARLIN lineage barcode sequence as complete when `cutadapt=True`.
4. Write barcode-matched FASTQs.

Important preprocessing parameters:

- `cutadapt`: whether to trim reads before extraction.
- `amp_cores`: parallelism for cutadapt and barcode extraction.
- `base_quality`: base-quality threshold.
- `linker1`, `linker2`, `mm_rate`: linker matching.
- `gzip_output`, `gzip_after_preprocess`: output compression behavior.

### 5.2 DARLIN Correction

After preprocessing, the shell script calls:

```text
script/Quality_Control/python/amplicon.py
```

The correction workflow:

1. Collapse raw molecules by barcode, UMI, and lineage barcode.
2. Filter low-read raw molecules with `--initial_reads_cutoff`.
3. Correct spatial barcodes to the whitelist.
4. Correct UMIs within each spatial barcode using `--umi_hd_threshold`.
5. Correct lineage barcodes using `--lb_error_rate` and `--lb_min_hd`.
6. Keep the major LR per SR/UR group using `--major_fraction_threshold_molecule`.
7. Compute SR-level reads-per-UMI slope `k = n_reads / n_UR`.
8. Filter low-quality SR groups with `--slope_cutoff`.
9. Filter final rows with `--final-reads-cutoff`.

Important columns:

- `SR`: spatial barcode / spot identity
- `UR`: UMI after correction
- `LR`: corrected lineage barcode
- `reads`: read support
- `reads_fraction`: LR fraction within an SR/UR group
- `k`: SR-level reads-per-UMI slope
- `n_LR`: number of unique LR values within an SR

Important outputs per locus:

```text
amplicon/results/<sample_name>/<CA|RA|TA>/
├── final.csv
├── dbit.log
├── Reads_counts_heatmap.png
├── UMI_counts_heatmap.png
├── lr_per_sr_hist.png
├── reads_cutoff_qc.png
├── reads_fraction_qc.png
└── sr_reads_vs_umis.png
```

`final.csv` is the main per-locus clone-call table used by downstream tissue-filtered plotting.

## 6. Tissue-Filtered Visualization

Entry point:

```text
script/Quality_Control/plot.sh
```

This step combines image-derived tissue membership and retained cell-count
information with mRNA and/or amplicon spatial results. Spots are filtered by
`in_tissue`; cell counts remain available for cell-based summaries.

The primary input is passed on the command line; additional plotting paths are
set in the shared QC config:

- `cell_number_file`: appended to the dataset config by the image step; usually `image/filtered_results.csv`
- `tissue_mask_file`: whole-image binary mask, usually `image/tissue_mask.png`
- `mrna_dir`: STARsolo `GeneFull` directory
- `amp_dir`: amplicon result directory
- `gray_path`: grayscale image for merged overlays
- the barcode whitelist is selected automatically by the chip argument

For mRNA data, the script calls:

```text
script/Quality_Control/python/mrna_cell.py
```

Key mRNA output:

```text
<mrna_dir>/raw/data_tissuefiltered.csv
<mrna_dir>/raw/umap_filtered.png
<mrna_dir>/raw/umi_filtered.png
<mrna_dir>/raw/gene_filtered.png
<mrna_dir>/raw/merged_umap_filtered.png
<mrna_dir>/filtered/data_tissuefiltered.csv
<mrna_dir>/filtered/umap_filtered.png
<mrna_dir>/filtered/umi_filtered.png
<mrna_dir>/filtered/gene_filtered.png
<mrna_dir>/filtered/merged_umap_filtered.png
```

For amplicon data, the script calls:

```text
script/Quality_Control/python/amplicon_cell.py
```

Key amplicon outputs are written per locus:

```text
<amp_dir>/<CA|RA|TA>/
├── tissuefiltered.csv
├── umi_filtered.png
└── merged_umi_filtered.png
```

When `gray.png` is available, `merge_on_gray.py` transforms each `*_filtered.png` image according to `--orientation` and `--swap_xy`, resizes the grayscale background, and composites the overlay on top.

## 7. Clone Analysis

Entry point:

```text
script/dbit.sh clone
script/Clone_Analysis/top_lr_pipeline.sh
```

The clone-analysis workflow reads its parameters from the config file set by
earlier steps. It expects the amplicon tissue-filtered outputs (from `amp_dir`):

```text
<input_dir>/
├── CA/tissuefiltered.csv
├── RA/tissuefiltered.csv
└── TA/tissuefiltered.csv
```

It also requires (from `bank_dir` and `cluster_csv` in config):

- an allele-bank directory containing `allele_bank_Gr_CA.csv.gz`, `allele_bank_Gr_RA.csv.gz`, and `allele_bank_Gr_TA.csv.gz`
- an mRNA cluster CSV with `x`, `y`, `leiden`, and optionally `color`

### 7.1 Allele-Bank Filter

Implemented in:

```text
script/Clone_Analysis/python/allele_bank_filter.py
```

This step:

1. Reads each label's `tissuefiltered.csv`.
2. Runs `darlin_core.analyze_sequences` on unique LR sequences.
3. Compares the resulting mutation strings against the label-specific allele-bank file.
4. Keeps LR rows whose mutation patterns are not present in the allele bank.
5. Writes a cache file named `.analyzed_cache.csv` under each label output directory.

Label-to-bank mapping:

```text
CA -> allele_bank_Gr_CA.csv.gz, config Col1a1
RA -> allele_bank_Gr_RA.csv.gz, config Rosa
TA -> allele_bank_Gr_TA.csv.gz, config Tigre
```

Output per label:

```text
tissuefiltered.bank_filtered.csv
```

### 7.2 Top LR Spatial Plot

Implemented in:

```text
script/Clone_Analysis/python/top_lr_plot.py
```

The plotter:

1. Reads `tissuefiltered.bank_filtered.csv`.
2. Groups by LR and ranks clones by:
   - number of unique SR spots
   - number of unique UR values
   - total read support
   - LR sequence
3. Selects the top `clone_top_n` LR clones.
4. Draws a Leiden cluster background from the mRNA cluster CSV.
5. Overlays hollow circles on spots containing the selected LR.
6. Sizes circles by unique UR count:
   - `1`
   - `2`
   - `3+`

Plot details:

- `--cluster-alpha` controls both Leiden background opacity and the cluster legend opacity.
- `--orientation <mode>` and `--swap-xy` control spot coordinate orientation, matching the QC pipeline conventions.
- Titles are formatted into a fixed-height area so output PNG dimensions stay constant.
- Edge padding and unclipped LR circles are used so boundary circles are not cut off.

Outputs per label:

```text
top_lr_plots/
├── topLR_001_srXXX_urXXX.png
├── topLR_002_srXXX_urXXX.png
└── <LABEL>_top_lr_plot_manifest.csv
```

The manifest contains:

- `rank`
- `LR`
- `unique_SR`
- `unique_UR`
- `total_reads`
- `plot_file`

## 8. Output Summary

Typical high-level result layout:

```text
sample_name/
├── transcriptome/
│   └── results/
├── image/
│   ├── filtered_results.csv
│   ├── tissue_mask.png
│   └── gray.png
└── amplicon/
    └── results/
        └── sample_name/
            ├── CA/
            ├── RA/
            └── TA/
```

Clone-analysis output layout:

```text
<output_dir>/
├── CA/
│   ├── tissuefiltered.n_LR_gt_count.csv
│   ├── tissuefiltered.n_LR_le_count.csv
│   ├── tissuefiltered.count_summary.txt
│   ├── tissuefiltered.bank_filtered.csv
│   ├── .analyzed_cache.csv
│   └── top_lr_plots/
│       ├── topLR_001_srXXX_urXXX.png
│       └── CA_top_lr_plot_manifest.csv
├── RA/
└── TA/
```

## 9. Recommended Debug Checks

1. Inspect preprocessing logs when barcode output is unexpectedly small:

   ```text
   <output_path>/<sample_name>_preprocess.log
   ```

2. Inspect STAR logs when mRNA matrices are missing:

   ```text
   results/<sample_name>/STAR.log
   results/<sample_name>/Solo.out/qc.log
   ```

3. Inspect image registration before trusting tissue-filtered plots:

   ```text
   image/result.png
   image/filtered_results.csv
   ```

4. Inspect amplicon correction summaries and QC plots:

   ```text
   amplicon/results/<sample>/<label>/dbit.log
   amplicon/results/<sample>/<label>/reads_fraction_qc.png
   amplicon/results/<sample>/<label>/sr_reads_vs_umis.png
   ```

5. Inspect clone-analysis intermediate files before interpreting top LR plots:

   ```text
   tissuefiltered.count_summary.txt
   tissuefiltered.bank_filtered.csv
   top_lr_plots/*_top_lr_plot_manifest.csv
   ```
