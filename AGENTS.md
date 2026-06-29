# Repository Guidelines

## Project Overview

DBiT-spatial-DARLIN is a quality-control pipeline for DBiT spatial DARLIN data. The main user-facing entry point is `script/dbit.sh`; worker scripts and Python helpers live under `script/Quality_Control/`. Clone analysis utilities live under `script/Clone_Analysis/`.

Primary QC scripts:

- `script/dbit.sh`
- `script/Quality_Control/mrna.sh`
- `script/Quality_Control/image.sh`
- `script/Quality_Control/amplicon.sh`
- `script/Quality_Control/plot.sh`

Detailed user workflow is documented in `README.md` and `docs/TECHNICAL_DOCUMENTATION.md`.

## Environment

Use `pixi` for environment management.

- Prefer the current project's pixi environment at `pixi.toml` when running commands from this repository.
- The relevant pixi environments are `default` and `image`.
- Run command examples from this repository with the project manifest when possible, for example:
  - `pixi run --manifest-path pixi.toml -e default <command>`
  - `pixi run --manifest-path pixi.toml -e image <command>`
- Do not modify `pixi.lock`.
- Avoid running environment solve/update commands unless explicitly requested.

## Development Notes

- Keep changes scoped to the requested pipeline behavior.
- Preserve the shell-script entry points and existing directory layout.
- Prefer `rg` for searching files and text.
- Before editing, check `git status --short`; this repository may already contain user changes.
- Do not revert user changes or unrelated modified files.

## Validation

There is no dedicated test suite visible in this repository. Validate changes with the narrowest relevant checks available, such as:

- shell syntax checks: `bash -n script/Quality_Control/<script>.sh`
- Python syntax checks through pixi: `pixi run --manifest-path pixi.toml -e default python -m py_compile <file>`
- image-related Python checks through `image` when they depend on TensorFlow, StarDist, OpenCV, or image codecs.

When validating commands that need project dependencies, use the current project pixi manifest first and avoid any command that would rewrite `pixi.lock`.
