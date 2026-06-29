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
├── dbit.sh
├── config.sh
├── Quality_Control/
│   ├── mrna.sh
│   ├── image.sh
│   ├── amplicon.sh
│   ├── plot.sh
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

## 3. Recommended Workflow

Run commands from the repository root:

```bash
cd /path/to/DBiT-spatial-DARLIN/script
```

Copy `config.sh` into the dataset directory, or keep the copy in the current
`script/` directory under a sample-specific name such as `config.sample.sh`.
Then edit its reference settings and pipeline parameters, and select
`execution_mode=local` or `execution_mode=hpc`.

```bash
# Option 1: store the config with the dataset
cp config.sh /path/to/sample_name/config.sh
bash dbit.sh mrna --input /path/to/transcriptome/fastq --config /path/to/sample_name/config.sh --chip 50-20

# Option 2: keep a renamed copy in the script directory
cp config.sh config.sample.sh
bash dbit.sh mrna --input /path/to/transcriptome/fastq --config config.sample.sh --chip 50-20
```

The mRNA step also accepts optional filtering thresholds. Supplied values are
appended to the dataset config; omitted values keep the config defaults:

```bash
bash dbit.sh mrna \
    --input /path/to/transcriptome/fastq \
    --config /path/to/config.sh \
    --chip 50-20 \
    --umi-min 900 \
    --gene-min 300 \
    --min-cell 3
```

These filtering options are accepted only by the `mrna` step.

The amplicon step provides its own optional filtering thresholds:

```bash
bash dbit.sh amplicon \
    --input /path/to/amplicon/fastq \
    --config /path/to/config.sh \
    --chip 50-20 \
    --initial-reads-cutoff 100 \
    --major-fraction-threshold-molecule 0.8 \
    --reads-cutoff 10 \
    --slope-cutoff 10
```

These options are accepted only by the `amplicon` step and are appended to the
dataset config when supplied.

Image requires both orientation parameters and stores them for reuse by plot:

```bash
bash dbit.sh image \
    --input /path/to/align.png \
    --config /path/to/config.sh \
    --chip 50-20 \
    --orientation horizontal \
    --swap-xy True
```

These required options are accepted only by the `image` step. The plot step
reads the values previously written to the dataset config and fails if they are
absent.

Plot is the final QC step. It reads all accumulated
paths and orientation settings from the dataset config, so it needs no input:

```bash
bash dbit.sh plot --config /path/to/config.sh
```

The recommended order is `mrna → amplicon → image → plot`. For `mrna`,
`amplicon`, and `image`, `--input` supports shell path completion and each run
appends its input and derived result paths to the dataset config.

`--chip` selects one of three presets defined in `dbit.sh`:

- `50-50`: 50x50 spots, spot length 50, interval 50
- `50-20`: 50x50 spots, spot length 20, interval 20
- `100-20`: 100x100 spots, spot length 20, interval 20

The selected chip name is appended to the per-dataset config. After it has been
written once, later steps may omit `--chip`; supplying it again appends a newer
selection. Grid dimensions remain defined only in `dbit.sh`.

The matching 50- or 100-barcode whitelist is selected automatically. Invoke
`dbit.sh` once per step that should be run or submitted.

The unified launcher passes the resolved chip settings to each independent QC
job. Chip presets are not duplicated in the config or worker scripts.

### Clone Analysis

```bash
bash Clone_Analysis/top_lr_pipeline.sh \
    -i /path/to/sample_name/amplicon/results/sample_name \
    -b /path/to/allele_bank \
    --cluster-csv /path/to/sample_name/transcriptome/results/sample_name/Solo.out/GeneFull/raw/data_cellfiltered.csv
```

This workflow:

1. splits clone calls by whether `n_LR > predicted cell count`;
2. filters LR sequences against label-specific allele-bank files;
3. plots the top LR clones on the mRNA Leiden cluster background.

## 4. Orientation

`image.sh` and `plot.sh` read the same orientation values from
the QC config:

```bash
orientation=normal       # normal, horizontal, vertical, or rotate
swap_xy=False            # True swaps the x and y axes
```

Common meanings:

- `normal`: no coordinate transform
- `horizontal`: left-right flip
- `vertical`: top-bottom flip
- `rotate`: 180-degree rotation
- `orientation=horizontal` with `swap_xy=True`: 90 degrees counterclockwise
- `orientation=vertical` with `swap_xy=True`: 90 degrees clockwise

Use the same orientation settings for image splitting and filtered-plot merging whenever possible. See [ORIENTATION.md](docs/ORIENTATION.md) for schematic examples.

`top_lr_plot.py` has its own `--rotate <0|90|180|270>` option because the clone-analysis plot orientation is handled separately from image merging.

## 5. Script Interfaces

Use `-h` or `--help` for the full current interface:

```bash
bash dbit.sh -h
bash dbit.sh mrna -h
bash dbit.sh amplicon -h
bash dbit.sh image -h
bash dbit.sh plot -h
bash Quality_Control/mrna.sh -h
bash Quality_Control/image.sh -h
bash Quality_Control/amplicon.sh -h
bash Quality_Control/plot.sh -h
bash Clone_Analysis/top_lr_pipeline.sh -h
```

Use `dbit.sh <step> --input <path> --config <file> [--chip <name>]` for data
steps and `dbit.sh plot --config <file> [--chip <name>]` for plotting.
Use `dbit.sh <step> -h` for step-specific command-line parameters; base
pipeline and SLURM settings are documented inline in the config template.

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
