# DBiT-spatial-DARLIN

Version 0.2.0

DBiT-spatial-DARLIN is a quality-control pipeline for DBiT spatial
transcriptome, registered tissue image, and DARLIN data.

For implementation details, see the
[technical documentation](docs/TECHNICAL_DOCUMENTATION.md).

## Expected data organization

Organize each dataset as follows:

```text
sample_name/
├── dbit.config.sh
├── transcriptome/
│   └── fastq/
│       ├── <sample>_<lane1>_R1.fq.gz
│       ├── <sample>_<lane1>_R2.fq.gz
│       ├── <sample>_<lane2>_R1.fq.gz
│       └── <sample>_<lane2>_R2.fq.gz
├── image/
│   ├── <sample>.jpg
│   └── mask.png
└── darlin/
    └── fastq/
        ├── <sample>_CA_R1.fq.gz
        ├── <sample>_CA_R2.fq.gz
        ├── <sample>_RA_R1.fq.gz
        ├── <sample>_RA_R2.fq.gz
        ├── <sample>_TA_R1.fq.gz
        └── <sample>_TA_R2.fq.gz
```

The DARLIN filenames must contain `CA`, `RA`, or `TA` so the locus can be
identified. The mRNA FASTQ directory may contain one or more `*_R1.fq.gz` and
`*_R2.fq.gz` pairs, but all pairs in that directory must be lanes or chunks
from the same biological library. They are combined into one STARsolo result.

The full-resolution image is passed to the image step. Its adjacent `mask.png`
must have the same canvas dimensions; the mask's nonzero/opaque region defines the DBiT frame.

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

Initialize the configuration for each dataset:

```bash
cd /path/to/sample_name
dbit init
```

This copies a configuration file to the current directory as `dbit.config.sh`.
The repository provides two templates in `config/`:

- `dbit.config.sh` — pre-filled configuration with common default values (used
  by `dbit init` when available)
- `dbit.config.example.sh` — blank template with all fields commented out

Edit the copied configuration, including `genome_dir`, execution mode, and
SLURM resources where applicable. Use `execution_mode=local` for a local run or
`execution_mode=hpc` for SLURM submission.

Run `dbit` from the dataset directory. By default it loads `./dbit.config.sh`;
use `--config <file>` only when the configuration is stored elsewhere.

## Usage

The recommended order is:

```text
mrna → saturation → darlin → image
```

Input paths and chip selections are stored in the dataset config and reused by
later commands. Current options and filtering thresholds are listed by the
command-line help:

```bash
dbit -h
dbit mrna -h
```

## Output organization

After completing the workflow, the main outputs are organized as follows:

```text
sample_name/
├── dbit.config.sh
├── transcriptome/
│   ├── fastq/
│   ├── saturation/
│   ├── fastq_umi_barcode/
│   ├── matrix/
│   └── results/
├── image/
│   ├── fullres_grayscale.png
│   ├── tissue_mask.png
│   └── tissue_positions.tsv.gz
└── darlin/
    ├── fastq/
    ├── fastq_umi_barcode/
    └── results/
```

The dataset config is updated with the resolved input and result paths so later
steps can reuse outputs from earlier steps.
