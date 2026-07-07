# DBiT-spatial-DARLIN Technical Documentation

This document describes the implementation details of the DBiT-spatial-DARLIN processing workflow. For runnable examples and a shorter overview, see the repository [README](../README.md).

## 1. Pipeline Scope

The pipeline contains five user-facing steps:

- `mrna`: preprocess transcriptome FASTQs, run STARsolo, and perform spatial QC.
- `image`: split the registered image and count cells with StarDist.
- `amplicon`: process DARLIN amplicon FASTQs and generate clone-call tables.
- `filter`: apply tissue filtering and merge spatial plots with the image.
- `clone`: filter LR sequences against allele banks and plot the top clones.

The corresponding shell entry points are:

```text
script/dbit.sh
script/Quality_Control/mrna.sh
script/Quality_Control/image.sh
script/Quality_Control/amplicon.sh
script/Quality_Control/filter.sh
script/Clone_Analysis/clone.sh
```

`dbit.sh` launches one step locally or through SLURM. It stores resolved input
and result paths, chip selection, orientation, and selected command-line
overrides in a per-dataset config so later steps can reuse them. Chip grid
dimensions and whitelist selection are resolved centrally by the launcher and
exported to worker scripts. Python commands run through the appropriate Pixi
environment.

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

`image.sh`, `filter.sh`, and clone analysis share orientation controls:

The shared QC config sets `orientation` to `normal`, `horizontal`, `vertical`,
or `rotate`; `swap_xy=True` additionally swaps the coordinate axes.

These parameters are documented in detail in [ORIENTATION.md](ORIENTATION.md). The clone-analysis pipeline reads `orientation` and `swap_xy` from the same shared config, matching the QC pipeline conventions.

Clone `--rotate` is separate from orientation alignment. It rotates only the
final top-LR spatial grid clockwise by `0`, `90`, `180`, or `270` degrees for
presentation; it does not rotate the photographed image or the plot legends.

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
--soloCellFilter None
--soloFeatures GeneFull
```

STARsolo cell calling is disabled. Downstream mRNA QC uses only the `raw`
matrix and applies the configured UMI, gene-count, and minimum-spot filters
before clustering.

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
├── pearson_residuals_hvg.h5ad
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

`clustered.h5ad` keeps the full raw count matrix, spot metadata, spatial
coordinates, dimensional reductions, and clustering results.
`pearson_residuals_hvg.h5ad` stores the Pearson-residual-normalized
highly-variable-gene matrix in `X`, with matching spot and gene metadata. The
normalized matrix is not duplicated inside `clustered.h5ad`.

`data_tissuefiltered.csv` is produced later by `filter.sh` after applying the image-derived tissue mask to mRNA spot data.

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

`filtered_results.csv` is the image-derived table used by `filter.sh`. It retains
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

Amplicon preprocessing and correction output is displayed in the terminal
while also being written to `<sample>_preprocess.log` and `dbit.log`.

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

`reads_fraction_mode` selects whether the major-LR fraction denominator uses
the sum of LR read counts (`sum`) or the maximum LR read count (`max`).

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
script/Quality_Control/filter.sh
```

This step combines image-derived tissue membership and retained cell-count
information with mRNA and/or amplicon spatial results. Spots are filtered by
`in_tissue`; cell counts remain available for cell-based summaries.

All filtering paths and the chip selection are read from the shared dataset
config; the `filter` step does not accept `--input` or `--chip`:

- `cell_number_file`: appended to the dataset config by the image step; usually `image/filtered_results.csv`
- `tissue_mask_file`: whole-image binary mask, usually `image/tissue_mask.png`
- `mrna_dir`: STARsolo `GeneFull` directory
- `amp_dir`: amplicon result directory
- `gray_path`: grayscale image for merged overlays
- the barcode whitelist is selected automatically from the stored chip

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

When `gray.png` is available, `filter.sh` passes the expected mRNA and amplicon
frame paths explicitly to `merge_on_gray.py`. The script transforms only those
named files according to `--orientation` and `--swap_xy`, resizes the grayscale
background, and composites the overlay on top. It does not search directories
for matching images.

The mRNA and amplicon filtered-plot commands write to `filtered_plot.log` and
print the same output to the terminal.

## 7. Clone Analysis

Entry point:

```text
script/Clone_Analysis/clone.sh
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

- an allele-bank directory containing one `.csv`, `.csv.gz`, `.tsv`, or
  `.tsv.gz` file for each locus; filenames must contain the corresponding
  uppercase `CA`, `RA`, or `TA` label
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
CA filename label -> config Col1a1
RA filename label -> config Rosa
TA filename label -> config Tigre
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
- Transformations are applied in this order: stored `orientation`, stored
  `swap_xy`, then clone-only `--rotate`.
- `--orientation <mode>` and `--swap-xy` control spot coordinate orientation, matching the QC pipeline conventions.
- `--rotate <0|90|180|270>` applies an additional clockwise rotation to the spatial grid for presentation without rotating the legends.
- `--x-spots-number` and `--y-spots-number` receive the complete chip grid
  dimensions resolved by `dbit.sh`. Clone and cluster coordinates are rejected
  if they fall outside that grid.
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
<amp_dir>/
├── CA/
│   ├── tissuefiltered.csv
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
   tissuefiltered.bank_filtered.csv
   top_lr_plots/*_top_lr_plot_manifest.csv
   ```
