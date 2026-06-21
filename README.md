## DBiT-spatial-DARLIN

Version 0.2.0

This repository contains scripts for DBiT spatial DARLIN data processing and quality control across three data types:

- transcriptome FASTQ data
- registered tissue images
- DARLIN amplicon FASTQ data

It also includes a clone-analysis workflow that filters DARLIN clone calls and plots top LR clones on an mRNA Leiden-cluster background.

Detailed references:

- [Technical documentation](docs/TECHNICAL_DOCUMENTATION.md)
- [Orientation handling](docs/ORIENTATION.md)

## 1. Project Layout

Recommended input layout:

```text
sample_name/
├── amplicon/
│   └── fastq/
│       ├── *_CA_R1.fq.gz
│       ├── *_CA_R2.fq.gz
│       ├── *_RA_R1.fq.gz
│       ├── *_RA_R2.fq.gz
│       ├── *_TA_R1.fq.gz
│       └── *_TA_R2.fq.gz
├── image/
│   ├── align.png
│   └── gray.png
└── transcriptome/
    └── fastq/
        ├── *_R1.fq.gz
        └── *_R2.fq.gz
```

Main script directories:

```text
script/
├── Quality_Control/
│   ├── dbit_mrna.sh
│   ├── image.sh
│   ├── dbit_amplicon.sh
│   ├── plot_cell_filtered.sh
│   └── python/
└── Clone_Analysis/
    ├── top_lr_pipeline.sh
    └── python/
```

## 2. Environment

This project uses `pixi`.

```bash
cd /path/to/DBiT-spatial-DARLIN
pixi install
```

The shell entry points call `pixi run -e ...` internally. You do not need to activate the environment manually.

Pixi environments:

- `default`: transcriptome, amplicon, filtered plotting, and clone analysis
- `image`: image splitting and StarDist segmentation

Useful environment options accepted by the shell scripts:

```bash
--pixi_env <name>
--pixi_env_dir /path/to/DBiT-spatial-DARLIN
```

## 3. Recommended Workflow

Run commands from the repository root:

```bash
cd /path/to/DBiT-spatial-DARLIN
```

### Step 1: Transcriptome QC

```bash
bash script/Quality_Control/dbit_mrna.sh \
    -f /path/to/sample_name/transcriptome/fastq \
    -w docs/barcodes/barcodes.tsv \
    --genome_dir /path/to/STAR_genome_index
```

This processes transcriptome FASTQs, runs STAR/Solo, and generates mRNA QC outputs.

### Step 2: Image Segmentation

Register the transcriptome plot to the tissue image first, crop the corresponding region, and save it as `align.png`.

```bash
bash script/Quality_Control/image.sh \
    -i /path/to/sample_name/image/align.png \
    --orientation normal
```

This splits the image by the DBiT grid, runs StarDist, and writes `filtered_results.csv`.

### Step 3: Amplicon QC

```bash
bash script/Quality_Control/dbit_amplicon.sh \
    -f /path/to/sample_name/amplicon/fastq \
    -w docs/barcodes/barcodes.tsv \
    -c True
```

This processes DARLIN amplicon FASTQs and writes per-label clone-call tables.

### Step 4: Cell-Filtered Plots

```bash
bash script/Quality_Control/plot_cell_filtered.sh \
    -c /path/to/sample_name/image/filtered_results.csv \
    -m /path/to/sample_name/transcriptome/results/sample_name/Solo.out/GeneFull \
    -a /path/to/sample_name/amplicon/results/sample_name \
    -w docs/barcodes/barcodes.tsv \
    -g /path/to/sample_name/image/gray.png \
    --orientation normal
```

This generates cell-filtered transcriptome and/or amplicon plots. If `gray.png` is provided or found next to `filtered_results.csv`, merged overlay images are also produced.

### Step 5: Clone Analysis

```bash
bash script/Clone_Analysis/top_lr_pipeline.sh \
    -i /path/to/sample_name/amplicon/results/sample_name \
    -b /path/to/allele_bank \
    --cluster-csv /path/to/sample_name/transcriptome/results/sample_name/Solo.out/GeneFull/raw/data_cellfiltered.csv
```

This workflow:

1. splits clone calls by whether `n_LR > predicted cell count`;
2. filters LR sequences against label-specific allele-bank files;
3. plots the top LR clones on the mRNA Leiden cluster background.

## 4. Orientation

`image.sh` and `plot_cell_filtered.sh` share the same orientation interface:

```bash
--orientation normal
--orientation horizontal
--orientation vertical
--orientation rotate
--swap_xy
```

Common meanings:

- `normal`: no coordinate transform
- `horizontal`: left-right flip
- `vertical`: top-bottom flip
- `rotate`: 180-degree rotation
- `horizontal --swap_xy`: 90 degrees counterclockwise
- `vertical --swap_xy`: 90 degrees clockwise

Use the same orientation settings for image splitting and filtered-plot merging whenever possible. See [ORIENTATION.md](docs/ORIENTATION.md) for schematic examples.

`top_lr_plot.py` has its own `--rotate <0|90|180|270>` option because the clone-analysis plot orientation is handled separately from image merging.

## 5. Script Interfaces

Use `-h` or `--help` for the full current interface:

```bash
bash script/Quality_Control/dbit_mrna.sh -h
bash script/Quality_Control/image.sh -h
bash script/Quality_Control/dbit_amplicon.sh -h
bash script/Quality_Control/plot_cell_filtered.sh -h
bash script/Clone_Analysis/top_lr_pipeline.sh -h
```

### `dbit_mrna.sh`

Required:

- `-f, --fastq_path <dir>`: directory containing transcriptome `*_R1.fq.gz` and `*_R2.fq.gz`

Common options:

- `-o, --output_path <dir>`: output directory; default is the parent directory of `fastq_path`
- `-w, --whitelist <path>`: barcode whitelist
- `--genome_dir <path>`: STAR genome index
- `--umi_min <num>`: minimum UMI count per spot; default `900`
- `--gene_min <num>`: minimum gene count per spot; default `300`
- `--x_spots_number <num>` and `--y_spots_number <num>`: grid dimensions; default `50` and `50`

### `image.sh`

Required:

- `-i, --image_path <path>`: registered input image, usually `align.png`

Common options:

- `-r, --result_path <path>`: output directory; default is the input image directory
- `--orientation <mode>` and `--swap_xy`: grid orientation controls
- `--model <name>`: StarDist model; default `2D_versatile_fluo`
- `--prob_thresh <num>` and `--nms_thresh <num>`: StarDist detection thresholds
- `--cutoff <num>`: cell-count cutoff used by downstream filtering; default `100`

### `dbit_amplicon.sh`

Required:

- `-f, --fastq_path <dir>`: directory containing amplicon `*_R1.fq.gz` and `*_R2.fq.gz`

Common options:

- `-o, --output_path <dir>`: output directory; default is the parent directory of `fastq_path`
- `-w, --whitelist <path>`: barcode whitelist
- `-c, --cutadapt <bool>`: whether to run cutadapt trimming
- `--umi_hd_threshold <num>`: UMI correction threshold; default `1`
- `--initial_reads_cutoff <num>`: minimum reads per raw LB/SB/UB molecule; default `100`
- `--major_fraction_threshold <float>`: major LR fraction threshold; default `0.8`
- `--reads_cutoff <num>`: minimum reads per final row; default `10`

### `plot_cell_filtered.sh`

Required:

- `-c, --cell_number_file <path>`: cell number file, usually `filtered_results.csv`
- at least one of `-m/--mrna_dir` or `-a/--amp_dir`

Common options:

- `-m, --mrna_dir <path>`: mRNA result directory
- `-a, --amp_dir <path>`: amplicon result directory
- `-w, --whitelist <path>`: required when `--amp_dir` is provided
- `-g, --gray_path <path>`: grayscale image for merge; default is `gray.png` next to `filtered_results.csv`
- `--orientation <mode>` and `--swap_xy`: transform filtered images before merging

### `top_lr_pipeline.sh`

Required:

- `-i, --input-dir <dir>`: directory containing `CA`, `RA`, and/or `TA` subdirectories with `cellfiltered.csv`
- `-b, --bank-dir <dir>`: directory containing `allele_bank_Gr_CA.csv.gz`, `allele_bank_Gr_RA.csv.gz`, and `allele_bank_Gr_TA.csv.gz`
- `--cluster-csv <path>`: mRNA cluster CSV with `x`, `y`, `leiden`, and optionally `color`

Common options:

- `--labels <label...>`: labels to process; default `CA RA TA`
- `--output-dir <dir>`: output directory; default is `input-dir`
- `--min-sequence-length <num>`: minimum sequence length passed to DARLIN allele analysis; default `20`
- `--top-n <num>`: number of top LR plots per label; default `10`
- `--rotate <0|90|180|270>`: rotate top LR plots; default `0`

## 6. Outputs

Typical output layout:

```text
sample_name/
├── transcriptome/
│   └── results/
├── image/
│   ├── filtered_results.csv
│   └── gray.png
└── amplicon/
    └── results/
        └── sample_name/
            ├── CA/
            ├── RA/
            └── TA/
```

Clone-analysis outputs are written under the selected `--output-dir`:

```text
<output_dir>/
├── CA/
│   ├── cellfiltered.n_LR_gt_count.csv
│   ├── cellfiltered.n_LR_le_count.csv
│   ├── cellfiltered.count_summary.txt
│   ├── cellfiltered.bank_filtered.csv
│   └── top_lr_plots/
│       ├── topLR_001_srXXX_urXXX.png
│       └── CA_top_lr_plot_manifest.csv
├── RA/
└── TA/
```
