# DBiT-spatial-DARLIN

Version 0.2.0

DBiT-spatial-DARLIN is a quality-control and clone-analysis pipeline for DBiT
spatial transcriptome, registered tissue image, and DARLIN amplicon data.

For implementation details, see the
[technical documentation](docs/TECHNICAL_DOCUMENTATION.md). For coordinate
orientation examples, see [orientation handling](docs/ORIENTATION.md).

## Expected data organization

Organize each dataset as follows:

```text
sample_name/
├── config.sh
├── transcriptome/
│   └── fastq/
│       ├── <sample>_R1.fq.gz
│       └── <sample>_R2.fq.gz
├── image/
│   ├── align.png
│   └── gray.png
└── amplicon/
    └── fastq/
        ├── <sample>_CA_R1.fq.gz
        ├── <sample>_CA_R2.fq.gz
        ├── <sample>_RA_R1.fq.gz
        ├── <sample>_RA_R2.fq.gz
        ├── <sample>_TA_R1.fq.gz
        └── <sample>_TA_R2.fq.gz
```

The amplicon filenames must contain `CA`, `RA`, or `TA` so the locus can be
identified. Clone analysis also requires an allele-bank directory containing
one `.csv`, `.csv.gz`, `.tsv`, or `.tsv.gz` file per locus. Each bank filename
must contain the corresponding uppercase `CA`, `RA`, or `TA` label.

`align.png` is the cropped registered image produced by manually aligning the
tissue image with the spatial transcriptome result. It is used as the input to
the image QC step. `gray.png` is the corresponding grayscale tissue image used
as the background for merged spatial plots.

## Installation and configuration

The project uses [Pixi](https://pixi.prefix.dev/latest/installation/) for
environment and dependency management. Install Pixi first on Linux or macOS:

```bash
curl -fsSL https://pixi.sh/install.sh | sh
source ~/.bashrc
```

After installation, enter the repository and install the locked project
environments and the user-level `dbit` command:

```bash
cd /path/to/DBiT-spatial-DARLIN
pixi run init
source ~/.bashrc
```

Copy the configuration template for each dataset:

```bash
cp /path/to/DBiT-spatial-DARLIN/script/config.sh /path/to/sample_name/config.sh
```

Edit the copied configuration, including `genome_dir`, `bank_dir`, execution
mode, and SLURM resources where applicable. Use `execution_mode=local` for a
local run or `execution_mode=hpc` for SLURM submission.

## Usage

The recommended order is:

```text
mrna → image → amplicon → plot → clone
```

### mRNA QC

Preprocess transcriptome FASTQs, run STAR, calculate QC metrics, and generate
spatial expression and clustering results.

```bash
dbit mrna \
    --config /path/to/sample_name/config.sh \
    --input /path/to/sample_name/transcriptome/fastq \
    --chip 50-20
```

### Image QC

Segment the registered image, count cells per spatial spot, and generate the
tissue mask used by later filtering.

```bash
dbit image \
    --config /path/to/sample_name/config.sh \
    --input /path/to/sample_name/image/align.png \
    --orientation normal \
    --swap-xy False
```

### Amplicon QC

Preprocess DARLIN amplicon FASTQs, correct barcodes, filter lineage records,
and generate per-locus results.

```bash
dbit amplicon \
    --config /path/to/sample_name/config.sh \
    --input /path/to/sample_name/amplicon/fastq
```

### Tissue-filtered plots

Apply image-derived cell/tissue filtering to the mRNA and amplicon results and
merge spatial overlays with the grayscale tissue image.

```bash
dbit plot --config /path/to/sample_name/config.sh
```

### Clone analysis

Filter clone calls against the locus-specific allele banks and plot the top LR
clones over the mRNA Leiden-cluster background.

```bash
dbit clone --config /path/to/sample_name/config.sh
```

For presentation, rotate only the clone grid while leaving its legends fixed:

```bash
dbit clone --config /path/to/sample_name/config.sh --rotate 90
```

The first input path and chip selection are stored in the dataset config and
reused by later commands. Current options and filtering thresholds are listed
by the command-line help:

```bash
dbit -h
dbit mrna -h
dbit image -h
dbit amplicon -h
dbit plot -h
dbit clone -h
```

## Output organization

After completing the workflow, the main outputs are organized as follows:

```text
sample_name/
├── config.sh
├── transcriptome/
│   ├── fastq/
│   ├── fastq_umi_barcode/
│   └── results/
│       └── <sample>/
│           └── Solo.out/
├── image/
│   ├── filtered_results.csv
│   ├── tissue_mask.png
│   └── result.png
└── amplicon/
    ├── fastq/
    ├── fastq_umi_barcode/
    └── results/
        └── <sample>/
            ├── CA/
            │   ├── final.csv
            │   ├── tissuefiltered.csv
            │   ├── tissuefiltered.bank_filtered.csv
            │   └── top_lr_plots/
            ├── RA/
            └── TA/
```

The dataset config is updated with the resolved input and result paths so later
steps can reuse outputs from earlier steps.
