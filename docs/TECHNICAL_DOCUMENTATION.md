# DBiT-spatial-DARLIN Technical Documentation

This document describes the implementation details of the DBiT-spatial-DARLIN processing workflow. For runnable examples and a shorter overview, see the repository [README](../README.md).

## 1. Pipeline Scope

The pipeline contains these user-facing steps:

- `mrna`: preprocess transcriptome FASTQs, run STARsolo, and perform spatial QC.
- `saturation`: downsample transcriptome FASTQs and run mRNA QC at each fraction.
- `darlin`: process DARLIN FASTQs and generate lineage-call tables.
- `image`: create a tissue mask, filter spatial results, and generate registered plots.

The corresponding shell entry points are:

```text
script/dbit.sh
script/step/mrna.sh
script/step/saturation.sh
script/step/darlin.sh
script/step/image.sh
```

`dbit.sh` launches one step locally or through SLURM. It stores resolved input
and result paths, chip selection, and selected command-line
overrides in a per-dataset config so later steps can reuse them. Chip grid
dimensions and barcode A/B whitelist paths are resolved centrally by the
launcher and exported separately to worker scripts. Python commands run through
the appropriate Pixi environment.

## 2. Shared Concepts

### Spatial Coordinates

Spatial spot coordinates are stored as integer `row` and `col` columns. These
replace the former `x` and `y` output columns:

- `row` (barcode A) increases from top to bottom.
- `col` (barcode B) increases from right to left.
- `(row=0, col=0)` is at the upper-right corner.

### Barcode Structure

Transcriptome and DARLIN preprocessing both extract a 16 bp spatial barcode and a 10 bp UMI. The 16 bp spatial barcode is built from two 8 bp components.

![Barcode and UMI structure](image/barcode.png)

## 3. Transcriptome Workflow

Entry point:

```text
script/step/mrna.sh
```

### 3.1 Preprocessing

The input directory must contain exactly one `*_R1.fq.gz` file and its matching
`*_R2.fq.gz` file. Both `dbit.sh` and the mRNA worker validate this requirement
before processing. The worker then calls the shared preprocessing entry point:

```text
script/step/python/preprocess.py
```

Main operations:

1. Match linker sequences.
2. Extract barcode and UMI sequence.
3. Optionally correct barcode components against the whitelist.
4. Write barcode/UMI FASTQ files for downstream STARsolo processing.

Important parameters:

- `linker1`, `linker2`: expected linker sequences.
- `mm_rate`: mismatch rate for linker matching.
- Barcode correction accepts exact whitelist matches and unambiguous
  single-substitution matches (fixed Hamming distance of 1).
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
<output_path>/results/Solo.out/GeneFull/
```

When incomplete STAR outputs are cleaned before a rerun.

### 3.3 mRNA QC and Clustering

After STARsolo, `mrna.sh` calls:

```text
script/step/python/mrna.py
```

This script contains the mRNA-specific matrix QC, filtering, plotting, and
clustering logic. DARLIN processing is handled only by the DARLIN workflow.

The clustering workflow:

1. Load the STARsolo count matrix.
2. Save raw QC metrics before normalization.
3. Filter spots using UMI and gene thresholds.
4. Normalize with a Scanpy Pearson-residual workflow.
5. Run PCA.
6. Build an SNN graph from PCA coordinates.
7. Run Leiden clustering.
8. Write spatial cluster plots and tabular outputs.

Before Pearson-residual normalization, spots whose total count is zero across
the selected highly variable genes are removed. This prevents low-depth
fractions from producing NaN residuals and causing PCA to fail.

Key defaults:

- `--umi_min`: `900`
- `--gene_min`: `300`
- `--min_cells`: `3`
- Pearson residual `theta`: `100`
- Pearson residual `n_top_genes`: `3000`, capped by available genes
- PCA `n_comps`: `50`, capped by available spots and genes
- SNN `n_neighbors`: `30`
- SNN `n_pcs`: `20`, capped by available PCs
- Leiden `resolution`: `0.2`
- Leiden `random_state`: `42`

Important mRNA outputs:

```text
Solo.out/GeneFull/raw/
├── spatial_metrics.csv
├── frame_umap.png
├── umap_legend.png
├── umi_filtered.png
└── gene_filtered.png
```

`spatial_metrics.csv` contains the spot-level information needed by saturation
analysis and tissue filtering:

- `row`, `col`
- raw QC metrics such as `umi_count` and `gene_count`
- `leiden`
- `color`


In `frame_umap.png`, coordinate `(row=0, col=0)` is at the upper-right corner.
`row` increases from top to bottom, and `col` increases from right to left.

### 3.4 Saturation analysis

`dbit saturation` reuses `mrna_fastq_path` stored by the mRNA step. The single
`script/step/saturation.sh` worker uses `seqtk sample` with the same seed
for both reads of every FASTQ pair, then invokes the complete mRNA worker once
for each fraction. The same exactly-one-pair validation is applied before
downsampling. Outputs follow this layout:

```text
<mRNA FASTQ parent>/saturation/<fraction>/
├── fastq/
│   ├── <sample>_<fraction>_R1.fq.gz
│   └── <sample>_<fraction>_R2.fq.gz
├── fastq_umi_barcode/
└── results/Solo.out/GeneFull/
```


## 4. DARLIN Workflow

Entry point:

```text
script/step/darlin.sh
```

The `dbit darlin` command processes CA, RA, and TA DARLIN FASTQ files. The
script infers the locus from sample names containing `CA`, `RA`, or `TA`, and
accepts either `sample-CA` or `sample_CA` naming.

### 4.1 Preprocessing

`darlin.sh` performs the locus-specific two-pass cutadapt trimming itself,
then passes the trimmed FASTQ pair to the shared `python/preprocess.py` barcode
and UMI extractor. The mRNA workflow calls the same extractor directly and
does not run cutadapt.

Main operations:

1. Optionally run cutadapt.
2. Extract spatial barcode and UMI.
3. Treat the DARLIN lineage barcode sequence as complete when `cutadapt=True`.
4. Write barcode-matched FASTQs.

Important preprocessing parameters:

- `cutadapt`: whether to trim reads before extraction.
- `darlin_cores`: parallelism for cutadapt and barcode extraction.
- `base_quality`: base-quality threshold.
- `linker1`, `linker2`, `mm_rate`: linker matching.
- `gzip_output`, `gzip_after_preprocess`: output compression behavior.

DARLIN preprocessing and correction output is displayed in the terminal
while also being written to `<sample>_preprocess.log` and `dbit.log`.

### 4.2 DARLIN Correction

After preprocessing, the shell script calls:

```text
script/step/python/darlin.py
```

The same script contains the DARLIN-specific FASTQ parsing, barcode/UMI/LR
correction, QC plotting, and `final.csv` output.
Barcode A and barcode B whitelist files are passed separately; SB sequences
are interpreted in B+A order, with A mapped to row and B mapped to column.

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
- `row_bc`, `col_bc`: corrected barcode A/B components
- `row`, `col`: zero-based spatial coordinates derived from the A/B whitelists

Important outputs per locus:

```text
darlin/results/<CA|RA|TA>/
├── final.csv
├── dbit.log
├── lineage_bc_length.png
├── lr_per_sr_hist.png
├── reads_cutoff_qc.png
├── reads_fraction_qc.png
└── sr_reads_vs_umis.png
```

`final.csv` is the main per-locus clone-call table used by downstream
tissue-filtered plotting. Spatial coordinates are written during DARLIN
processing, so the later filtering step does not need a barcode whitelist.

## 5. Image Workflow

Entry point:

```text
script/step/image.sh
```

### 5.1 Tissue Segmentation

The image workflow accepts a full-resolution image. An adjacent `mask.png`
with identical dimensions defines the spatial frame: its alpha channel is used
when it contains both transparent and opaque pixels, otherwise its nonzero
grayscale region is used.

Main steps:

1. Locate the frame bounding box from the adjacent mask.
2. Crop the corresponding region from the full-resolution image.
3. Generate a binary tissue mask within that region from image intensity,
   local signal density, morphology, and connected-region size.
4. Evaluate each configured DBiT spot against the tissue mask and retain spots
   whose tissue coverage is at least the internal threshold.
5. Write the retained barcode indices and their center positions in the
   coordinate system of the original full-resolution image.

Implementation file:

```text
script/step/python/image_segment.py
```

`image_segment.py` contains frame detection, tissue-mask generation, spatial
grid placement, and tissue-spot filtering. It does not segment or count cells.

Important outputs:

```text
image/
├── tissue_mask.png
└── tissue_positions.tsv.gz
```

`tissue_mask.png` is cropped to the frame bounding box.
`tissue_positions.tsv.gz` is the image-derived index used by the image step and
contains `barcode`, `array_row`, `array_col`, `pxl_row_in_fullres`, and
`pxl_col_in_fullres`. Barcodes are written in B+A sequence order. Pixel
coordinates locate each retained spot center in the original image. Spots with
less than the internal minimum tissue coverage are omitted from this file.

Tissue-mask thresholds and cleanup parameters are internal defaults; no
additional user configuration is required.

### 5.2 Tissue Filtering and Visualization

After segmentation, the image step joins the retained image-derived spot index
with mRNA and/or DARLIN spatial results. Spots absent from
`tissue_positions.tsv.gz` are removed.

After generating the tissue-position index, the image worker filters any mRNA
or DARLIN result paths available in the shared dataset config:

- `tissue_positions_file`: retained spot index, usually `image/tissue_positions.tsv.gz`
- `tissue_mask_file`: frame-sized binary mask, usually `image/tissue_mask.png`
- `mrna_dir`: STARsolo `GeneFull` directory
- `darlin_dir`: DARLIN result directory
- `fullres_image_path`: original full-resolution image
- `frame_mask_file`: adjacent mask that selects the frame

For configured mRNA and/or DARLIN data, the image worker calls one filter
entry point:

```text
script/step/python/image_filter.py
```

Key mRNA output:

```text
<mrna_dir>/raw/spatial_metrics_tissuefiltered.csv
<mrna_dir>/raw/umap_filtered.png
<mrna_dir>/raw/umi_filtered.png
<mrna_dir>/raw/gene_filtered.png
<mrna_dir>/raw/merged_umap_filtered.png
```

The image step reads the uncompressed `matrix.mtx`, `barcodes.tsv`, and
`features.tsv` files from `GeneFull/raw`, then writes the tissue-filtered raw
count matrix in compressed 10x format while preserving the original feature
and barcode order:

```text
mrna/matrix/
├── matrix.mtx.gz
├── barcodes.tsv.gz
├── features.tsv.gz
├── tissue_positions.tsv.gz
└── <original-image-name>
```

The tissue-position file and original full-resolution image are copied into the
same directory so the filtered matrix and its spatial context can be moved as
one unit.

The only DARLIN filter output is written per locus:

```text
<darlin_dir>/<CA|RA|TA>/
└── tissuefiltered.csv
```

The image worker passes the expected mRNA frame paths explicitly to
`merge_on_image.py`.
The script crops the original image to the mask bounding box, resizes that crop
to the spatial-frame dimensions, and composites the overlay on top. It does not
search directories for matching images.

The mRNA filtered-plot command writes to `filtered_plot.log` and prints the
same output to the terminal.

## 6. Output Summary

Typical high-level result layout:

```text
sample_name/
├── transcriptome/
│   └── results/
├── darlin/
│   └── results/
│       ├── CA/
│       ├── RA/
│       └── TA/
└── image/
    ├── tissue_mask.png
    └── tissue_positions.tsv.gz
```

## 7. Recommended Debug Checks

1. Inspect preprocessing logs when barcode output is unexpectedly small:

   ```text
   <output_path>/<sample_name>_preprocess.log
   ```

2. Inspect STAR logs when mRNA matrices are missing:

   ```text
   results/STAR.log
   results/Solo.out/qc.log
   ```

3. Inspect DARLIN correction summaries and QC plots:

   ```text
   darlin/results/<label>/dbit.log
   darlin/results/<label>/reads_fraction_qc.png
   darlin/results/<label>/sr_reads_vs_umis.png
   ```

4. Inspect image registration before trusting tissue-filtered plots:

   ```text
   image/tissue_mask.png
   image/tissue_positions.tsv.gz
   ```
