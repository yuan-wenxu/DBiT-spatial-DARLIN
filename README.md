## DBiT-spatial-DARLIN

Version 0.2.0

This version includes quality control for DBiT spatial DARLIN data from three modalities: transcriptome, image, and amplicon.

Expected input layout:

```text
sample_name/
├── amplicon/
|   └── fastq/
|       ├── *_CA_R1.fq.gz
|       ├── *_CA_R2.fq.gz
|       ├── *_RA_R1.fq.gz
|       ├── *_RA_R2.fq.gz
|       ├── *_TA_R1.fq.gz
|       └── *_TA_R2.fq.gz
├── image/
|   └── your_image.tif
└── transcriptome/
    └── fastq/
        ├── *_R1.fq.gz
        └── *_R2.fq.gz
```

Script directory: `script`

Detailed technical documentation: [TECHNICAL_DOCUMENTATION.md](docs/TECHNICAL_DOCUMENTATION.md)

Orientation handling documentation: [ORIENTATION.md](docs/ORIENTATION.md)

Quality control entry points:

1. `script/Quality_Control/dbit_mrna.sh`: process transcriptome data.
2. `script/Quality_Control/image.sh`: split the registered image, run StarDist segmentation, and predict cell numbers.
3. `script/Quality_Control/dbit_amplicon.sh`: process amplicon data.
4. `script/Quality_Control/plot_cell_filtered.sh`: generate cell-filtered plots and optionally merge them onto `gray.png`.

## 1. Environment

This project uses `pixi` as the environment manager.

- `default` environment: used by `dbit_mrna.sh`, `dbit_amplicon.sh`, and `plot_cell_filtered.sh`
- `image` environment: used by `image.sh`

Environment files are included in the repository:

- `pixi.toml`
- `pixi.lock`

Setup:

```bash
cd /path/to/DBiT-spatial-DARLIN
pixi install
```

The shell scripts call `pixi run -e ...` internally, so environments do not need to be activated manually.

By default, each script looks for `pixi.toml` at the repository root. If commands are launched from another location or with a custom environment, use:

```bash
--pixi_env <name> --pixi_env_dir /path/to/DBiT-spatial-DARLIN
```

Validated with:

- `pixi 0.70.1`
- `pixi run -e default ...`
- `pixi run -e image ...`

Key tools and packages verified in the pixi environments:

- `default`: `STAR`, `samtools`, `seqtk`, `fastp`, `cutadapt`, `scanpy`, `umi_tools`
- `image`: `tensorflow`, `stardist`, `opencv-python-headless`, `imagecodecs`

## 2. Quick Start

Use `-h` or `--help` to inspect the current interface for each script:

```bash
bash script/Quality_Control/dbit_mrna.sh -h
bash script/Quality_Control/image.sh -h
bash script/Quality_Control/dbit_amplicon.sh -h
bash script/Quality_Control/plot_cell_filtered.sh -h
```

Typical commands:

```bash
# Transcriptome: output defaults to the parent directory of the fastq folder.
bash script/Quality_Control/dbit_mrna.sh \
    -f /path/to/sample_name/transcriptome/fastq \
    -w docs/barcodes/barcodes.tsv \
    --genome_dir /path/to/STAR_genome_index

# Image segmentation: result path is optional and defaults to the image folder.
bash script/Quality_Control/image.sh \
    -i /path/to/sample_name/image/align.png \
    -r /path/to/sample_name/image \
    --orientation normal

# Amplicon: output defaults to the parent directory of the fastq folder.
bash script/Quality_Control/dbit_amplicon.sh \
    -f /path/to/sample_name/amplicon/fastq \
    -w docs/barcodes/barcodes.tsv \
    -c True

# Cell-filtered plots: whitelist is required only when amplicon results are provided.
bash script/Quality_Control/plot_cell_filtered.sh \
    -c /path/to/sample_name/image/filtered_results.csv \
    -m /path/to/sample_name/transcriptome/results/sample_name/Solo.out/GeneFull \
    -a /path/to/sample_name/amplicon/results/sample_name \
    -w docs/barcodes/barcodes.tsv \
    -g /path/to/sample_name/image/gray.png \
    --orientation normal
```

Orientation options used by `image.sh` and `plot_cell_filtered.sh` are `normal`, `horizontal`, `vertical`, and `rotate`. Use `--swap_xy` together with `--orientation` for 90-degree rotation. See [ORIENTATION.md](docs/ORIENTATION.md) for examples and schematic images.

### Execution Order

1. Run `script/Quality_Control/dbit_mrna.sh` to process transcriptome FASTQ files. This generates transcriptome QC outputs, including plots that can be used for image registration.
2. Register the transcriptome plot to the ssDNA staining image, crop the corresponding region, save it as `align.png`, and generate a grayscale image named `gray.png`.
3. Run `script/Quality_Control/image.sh` on `align.png` to split the image, run StarDist segmentation, and generate `filtered_results.csv`.
4. Run `script/Quality_Control/dbit_amplicon.sh` to process amplicon FASTQ files. This can be completed before the final integration step.
5. Run `script/Quality_Control/plot_cell_filtered.sh` with `filtered_results.csv` to generate cell-filtered transcriptome and/or amplicon plots. If `gray.png` is available, merged images are generated automatically.

## 3. Script Interfaces

### `dbit_mrna.sh`

Required:

- `-f, --fastq_path <dir>`: directory containing `*_R1.fq.gz` and `*_R2.fq.gz`

Common options:

- `-o, --output_path <dir>`: output directory; defaults to the parent directory of `fastq_path`
- `-w, --whitelist <path>`: barcode whitelist file
- `--genome_dir <path>`: STAR genome index
- `--scratch <path>`: optional scratch directory for intermediate files
- `--pixi_env <name>`: pixi environment name; default `default`
- `--pixi_env_dir <path>`: directory containing `pixi.toml`; default repository root

Selected QC/spatial options:

- `--umi_min <num>`: minimum UMI count per spot; default `900`
- `--gene_min <num>`: minimum gene count per spot; default `300`
- `--min_cells <num>`: minimum cells per gene; default `3`
- `--x_spots_number <num>` and `--y_spots_number <num>`: grid dimensions; default `50` and `50`
- `--length_spot <num>`, `--interval <num>`, `--pixel_length <float>`: spot geometry parameters

### `image.sh`

Required:

- `-i, --image_path <path>`: input image, usually `align.png`

Common options:

- `-r, --result_path <path>`: output directory; defaults to the image file directory
- `--orientation <mode>`: `normal`, `horizontal`, `vertical`, or `rotate`; default `normal`
- `--swap_xy`: swap x/y axes after applying `--orientation`
- `--scratch <path>`: optional scratch directory for intermediate image files
- `--pixi_env <name>`: pixi environment name; default `image`
- `--pixi_env_dir <path>`: directory containing `pixi.toml`; default repository root

Selected image/segmentation options:

- `--x_spots_number <num>` and `--y_spots_number <num>`: grid dimensions; default `50` and `50`
- `--length_spot <num>`, `--interval <num>`, `--pixel_length <float>`: spot geometry parameters
- `--put_text <bool>` and `--font_size <num>`: label display options for split preview
- `--top_value <num>` and `--number_of_top_values <num>`: image quality-control thresholds
- `--model <name>`, `--prob_thresh <num>`, `--nms_thresh <num>`: StarDist model and detection thresholds
- `--cutoff <num>`: cell-count cutoff used by `cell_filter.py`; default `100`

### `dbit_amplicon.sh`

Required:

- `-f, --fastq_path <dir>`: directory containing amplicon `*_R1.fq.gz` and `*_R2.fq.gz`

Common options:

- `-o, --output_path <dir>`: output directory; defaults to the parent directory of `fastq_path`
- `-w, --whitelist <path>`: barcode whitelist file
- `-c, --cutadapt <bool>`: whether to perform cutadapt trimming; default `False`
- `--scratch <path>`: optional scratch directory for intermediate files
- `--pixi_env <name>`: pixi environment name; default `default`
- `--pixi_env_dir <path>`: directory containing `pixi.toml`; default repository root

Selected correction options:

- `--umi_hd_threshold <num>`: UMI clustering edit-distance threshold; default `1`
- `--lb_error_rate <float>`: per-base error rate for LB threshold; default `0.02`
- `--major_fraction_threshold <float>`: major LR fraction threshold; default `0.8`
- `--reads_cutoff <num>`: minimum supported reads; default `10`
- `--slope_cutoff <num>`: minimum reads/UMIs slope per spot; default `10`

### `plot_cell_filtered.sh`

Required:

- `-c, --cell_number_file <path>`: cell number file, usually `filtered_results.csv`
- At least one of `-m/--mrna_dir` or `-a/--amp_dir`

Common options:

- `-m, --mrna_dir <path>`: mRNA results directory
- `-a, --amp_dir <path>`: amplicon results directory
- `-w, --whitelist <path>`: barcode whitelist file; required when `--amp_dir` is provided
- `-g, --gray_path <path>`: grayscale image for merge; defaults to `gray.png` next to `cell_number_file`
- `--orientation <mode>`: transform filtered images before merging; default `normal`
- `--swap_xy`: swap x/y axes after applying `--orientation`
- `--pixi_env <name>`: pixi environment name; default `default`
- `--pixi_env_dir <path>`: directory containing `pixi.toml`; default repository root

Spatial options passed to both mRNA and amplicon cell-filtered plotting:

- `--x_spots_number <num>` and `--y_spots_number <num>`: grid dimensions; default `50` and `50`
- `--length_spot <num>`, `--interval <num>`, `--pixel_length <float>`: spot geometry parameters

## 4. Output Directory Structure

Main output locations:

```text
sample_name/
├── amplicon/
|   └── results/
├── image/
|   ├── filtered_results.csv
|   └── gray.png
└── transcriptome/
    └── results/
```

Key result example: clustering results registered to brain slices.

<p align="center">
    <img src="docs/image/test.png">
<p>
