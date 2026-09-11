# Repository Guidelines

## Project and Environment

DBiT-spatial-DARLIN is a QC pipeline whose user-facing entry point is
`script/dbit.sh`. Workflow workers and their Python programs live in
`script/step/`.
Workflow details belong in `README.md` and `docs/TECHNICAL_DOCUMENTATION.md`.

- Use the repository `pixi.toml` and run all code, including OpenCV and
  image-codec work, in `default`.
- Do not edit `pixi.lock` or solve/update environments unless requested.
- Preserve the shell entry points and existing directory layout.

## Script Conventions

- Put Python `parse_args()` immediately after imports. Expose only `--` options,
  including `--help`; do not add short aliases.
- Keep worker-shell help text centralized in `script/dbit.sh`; worker shell
  scripts should not duplicate user-facing help blocks.
- Keep each Python worker self-contained: define its I/O, coordinate, and
  plotting helpers in the script that uses them.
- Pass barcode A and B explicitly when both axes depend on whitelist order.
  Concatenated DBiT barcodes use B+A sequence order, with A mapped to row and B
  mapped to column.
- Keep changes scoped. Check `git status --short` first and preserve unrelated
  user changes.
- Use `rg` for searches. Convert an explicit Windows path directly, for example
  `E:\data\file.h5ad` to `/mnt/e/data/file.h5ad`, without searching `/mnt`.

## Plot Style

Specify colors directly in each plotting script; do not add a shared
`plot_style.py` module.

- General data: blue `#0072B2`.
- Secondary metrics: orange `#E69F00`; fitted curves may use green `#009E73`.
- Thresholds and segmentation boundaries: vermilion `#D55E00`.
- Continuous heatmaps use `Reds`; spatial intensity frames use red `#D73027`.
- Leiden colors are written by mRNA outputs.
- Plot titles use size 12; axis labels and tick labels use size 10.
- Preserve semantic neutral colors such as grayscale image backgrounds,
  transparent masks, and black outlines.

## Validation

Use the narrowest relevant check; there is no dedicated test suite.

- Shell: `bash -n <script>`.
- Python: `pixi run --manifest-path pixi.toml -e default python -m py_compile <file>`.
- OpenCV or image codecs: validate through `-e default`.
- For material plotting changes, generate a small example under `/tmp` and
  remove it after inspection.
