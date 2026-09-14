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

The input directory must contain one or more matching `*_R1.fq.gz` and
`*_R2.fq.gz` pairs. All pairs in one directory are treated as lanes or chunks
of the same biological library and are combined into one STARsolo result. Both
`dbit.sh` and the mRNA worker require a matching R2 for every discovered R1.
Additional R2 files without a matching R1 are ignored.
The worker calls the shared preprocessing entry point once per pair:

```text
script/step/python/preprocess.py
```

Main operations:

1. Match linker sequences.
2. Extract barcode and UMI sequence.
3. Optionally correct barcode components against the whitelist.
4. Write barcode/UMI FASTQ files for downstream STARsolo processing.

Each pair keeps separate preprocessing outputs and a manifest containing the
input file metadata and relevant preprocessing parameters. This supports
per-pair restart while invalidating stale outputs when an input or parameter
changes.

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

After every input pair has been preprocessed, the worker supplies the
comma-separated R2 files as the cDNA input and the corresponding
comma-separated R1 files as the final barcode/UMI input to one STAR invocation.
The input ordering is deterministic. `results/star_input_manifest.tsv` records
the complete ordered input set and key STARsolo parameters; complete STAR
outputs are reused only when that manifest still matches.

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

Incomplete or outdated STAR outputs are cleaned before a rerun.

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

`dbit saturation` reuses `mrna_fastq_path` stored by the mRNA step. The
`script/step/saturation.sh` worker uses `seqtk sample` with the same seed for
both reads of every FASTQ pair, then invokes the complete mRNA worker once for
each fraction. All downsampled pairs for a fraction are combined into one
STARsolo result. Outputs follow this layout:

```text
<mRNA FASTQ parent>/saturation/<fraction>/
├── fastq/
│   ├── <chunk1>_<fraction>_R1.fq.gz
│   ├── <chunk1>_<fraction>_R2.fq.gz
│   ├── <chunk2>_<fraction>_R1.fq.gz
│   └── <chunk2>_<fraction>_R2.fq.gz
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
with identical canvas dimensions defines the spatial frame: its alpha channel
is used when it contains both transparent and opaque pixels, otherwise its
nonzero grayscale region is used. The visible mask may be smaller than the
configured spatial frame when the frame extends beyond an image boundary. In
that case, the touched boundary and configured grid dimensions determine the
off-image part of the frame. A short mask extent that does not touch exactly one
boundary is rejected because its alignment is ambiguous.

Main steps:

1. Locate the visible frame bounding box from the adjacent mask and reconstruct
   any portion clipped by an image boundary.
2. Convert the complete source image to grayscale and save it at the original
   pixel dimensions.
3. Crop the corresponding region from the full-resolution image.
4. Generate a binary tissue mask within that region from image intensity,
   local signal density, morphology, and connected-region size.
5. Evaluate each configured DBiT spot against the tissue mask and assign
   `in_tissue=1` to tissue spots or `in_tissue=0` to non-tissue spots.
6. Write every spot's barcode, tissue status, and center position in the
   coordinate system of the original full-resolution image.

Implementation file:

```text
script/step/python/image_segment.py
```

`image_segment.py` contains frame detection, full-resolution grayscale export,
tissue-mask generation, spatial grid placement, and tissue-status assignment.
It does not segment or count cells.

Important outputs:

```text
image/
├── fullres_grayscale.png
├── tissue_mask.png
└── tissue_positions.tsv.gz
```

`fullres_grayscale.png` is a single-channel rendering of the complete source
image and has exactly the same width and height as that source image.
`tissue_mask.png` is cropped to the visible frame bounding box. Spots partly
beyond the image are evaluated with the unavailable area treated as
non-tissue. A spot whose center is outside the full-resolution image is always
assigned `in_tissue=0`, regardless of its visible tissue fraction.
`tissue_positions.tsv.gz` is the image-derived index used by the image step and
contains `barcode`, `in_tissue`, `array_row`, `array_col`,
`pxl_row_in_fullres`, and `pxl_col_in_fullres`, in that order. Barcodes are
written in B+A sequence order. Pixel coordinates locate every spot center in
the original image. Tissue spots use `in_tissue=1`; non-tissue spots use
`in_tissue=0`.

Tissue-mask thresholds and cleanup parameters are internal defaults; no
additional user configuration is required.

### 5.2 Tissue Filtering and Visualization

After segmentation, the image step selects rows with `in_tissue=1` and joins
them with mRNA and/or DARLIN spatial results. The `in_tissue` column is retained
in the filtered CSV outputs.

After generating the tissue-position index, the image worker filters any mRNA
or DARLIN result paths available in the shared dataset config:

- `tissue_positions_file`: all-spot tissue index, usually `image/tissue_positions.tsv.gz`
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
`features.tsv` files from `GeneFull/raw`, then copies their complete contents
to compressed 10x files without filtering matrix columns. Tissue membership is
recorded in `tissue_positions.tsv.gz` instead:

```text
mrna/matrix/
├── matrix.mtx.gz
├── barcodes.tsv.gz
├── features.tsv.gz
├── tissue_positions.tsv.gz
└── fullres_grayscale.png
```

The tissue-position file and full-resolution grayscale image are copied into
the same directory so the matrix and its spatial context can be moved as one
unit. The original color image is not copied into `mrna/matrix`.

The only DARLIN filter output is written per locus:

```text
<darlin_dir>/<CA|RA|TA>/
└── tissuefiltered.csv
```

The image worker passes the expected mRNA frame paths explicitly to
`merge_on_image.py`.
The script crops the original image to the visible mask bounding box. If the
configured spatial frame crosses an image boundary, it first crops the overlay
to the corresponding visible subsection (for example, a frame crossing the top
and right boundaries retains its lower-left subsection). It then resizes the
background crop to that subsection and composites the overlay on top. It does
not search directories for matching images.

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
    ├── fullres_grayscale.png
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
   image/fullres_grayscale.png
   image/tissue_mask.png
   image/tissue_positions.tsv.gz
   ```
