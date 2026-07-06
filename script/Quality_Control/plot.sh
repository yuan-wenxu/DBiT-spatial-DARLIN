#!/bin/bash
set -o pipefail

show_help() {
    cat << EOF
Usage: $0 <config_file>

Plot tissue-filtered results for mRNA and/or amplicon data while retaining cell counts.

Arguments:
  config_file   Per-dataset QC configuration file


Examples:
  $0 config.sh
EOF
}

SCRIPT_DIR=${QC_SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)} || exit 1
REPO_DIR=${REPO_DIR:-$(cd "$SCRIPT_DIR/../.." && pwd)} || exit 1
PYTHON_DIR="$SCRIPT_DIR/python"
MERGE_SCRIPT="$PYTHON_DIR/image_process/merge_on_gray.py"

if [[ ${1:-} == -h || ${1:-} == --help ]]; then show_help; exit 0; fi
if [[ $# -ne 1 ]]; then show_help >&2; exit 1; fi
config_file=$1
if [[ ! -f "$config_file" ]]; then
    echo "Error: config file not found: $config_file" >&2; exit 1
fi

source "$config_file"
pixi_env=${pixi_env:-default}
pixi_env_dir=${pixi_env_dir:-$REPO_DIR}

if [[ -z ${cell_number_file:-} ]]; then
    echo "Error: cell_number_file must be set in the QC config." >&2; exit 1
fi
if [[ -z ${orientation:-} || -z ${swap_xy:-} ]]; then
    echo "Error: orientation and swap_xy must first be written by the image step." >&2; exit 1
fi
case "${swap_xy,,}" in
    true) swap_xy=True ;;
    false) swap_xy=False ;;
    *) echo "Error: swap_xy must be True or False; got '$swap_xy'." >&2; exit 1 ;;
esac
for variable in x_spots_number y_spots_number length_spot interval; do
    if [[ -z ${!variable:-} ]]; then
        echo "Run this script through dbit.sh so --chip is resolved." >&2
        exit 1
    fi
done

normalize_dir_path() {
    local path="$1"
    while [[ "$path" != "/" && "$path" == */ ]]; do
        path="${path%/}"
    done
    printf '%s\n' "$path"
}

merge_with_gray() {
    local frame
    local merge_args
    if [[ -f "$gray_path" ]]; then
        (
            cd "$pixi_env_dir" || exit 1
            merge_args=(
                --gray "$gray_path"
                --orientation "$orientation"
            )
            for frame in "$@"; do
                merge_args+=(--frame "$frame")
            done
            if [[ "$swap_xy" == True ]]; then
                merge_args+=(--swap_xy)
            fi
            pixi run -e "$pixi_env" python "$MERGE_SCRIPT" "${merge_args[@]}"
        )
    else
        echo "Warning: gray image not found at $gray_path, skipping merge."
    fi
}

# Validate inputs
if [ -z "$mrna_dir" ] && [ -z "$amp_dir" ]; then
    echo "Error: at least one of mrna_dir or amp_dir must be set in config" >&2
    exit 1
fi
if [ -n "$amp_dir" ] && [ -z "$whitelist_path" ]; then
    echo "Error: whitelist_path is required when amp_dir is set" >&2
    exit 1
fi

if [ -n "$mrna_dir" ]; then
    mrna_dir=$(normalize_dir_path "$mrna_dir")
fi
if [ -n "$amp_dir" ]; then
    amp_dir=$(normalize_dir_path "$amp_dir")
fi

if [ -z "$gray_path" ]; then
    gray_dir=$(normalize_dir_path "$(dirname "$cell_number_file")")
    gray_path="$gray_dir/gray.png"
fi

case "$orientation" in
    normal|horizontal|vertical|rotate)
        ;;
    *)
        echo "Error: orientation must be one of normal, horizontal, vertical, rotate" >&2
        exit 1
        ;;
esac

if [[ ! -d "$pixi_env_dir" ]]; then
    echo "Error: pixi environment dir does not exist: $pixi_env_dir" >&2
    exit 1
fi

# mRNA tissue-filtered plot (only when mrna_dir is provided)
if [ -n "$mrna_dir" ]; then
    (
        cd "$pixi_env_dir" || exit 1
        pixi run -e "$pixi_env" python "$PYTHON_DIR/mrna_cell.py" \
            -c "$cell_number_file" \
            -d "$mrna_dir" \
            --x_spots_number "$x_spots_number" \
            --y_spots_number "$y_spots_number" \
            --length_spot "$length_spot" \
            --interval "$interval" \
            --pixel_length "$pixel_length"
    ) 2>&1 | tee "$mrna_dir/filtered_plot.log" || exit 1
    merge_with_gray \
        "$mrna_dir/raw/umap_filtered.png" \
        "$mrna_dir/raw/umi_filtered.png" \
        "$mrna_dir/raw/umi_per_cell_filtered.png" \
        "$mrna_dir/raw/gene_filtered.png" \
        "$mrna_dir/raw/gene_per_cell_filtered.png" || exit 1
fi

# Amplicon tissue-filtered plot (only when amp_dir is provided)
if [ -n "$amp_dir" ]; then
    (
        cd "$pixi_env_dir" || exit 1
        pixi run -e "$pixi_env" python "$PYTHON_DIR/amplicon_cell.py" \
            -c "$cell_number_file" \
            -d "$amp_dir" \
            -w "$whitelist_path" \
            --x_spots_number "$x_spots_number" \
            --y_spots_number "$y_spots_number" \
            --length_spot "$length_spot" \
            --interval "$interval" \
            --pixel_length "$pixel_length"
    ) 2>&1 | tee "$amp_dir/filtered_plot.log" || exit 1
    merge_with_gray \
        "$amp_dir/CA/umi_filtered.png" \
        "$amp_dir/CA/umi_per_cell_filtered.png" \
        "$amp_dir/RA/umi_filtered.png" \
        "$amp_dir/RA/umi_per_cell_filtered.png" \
        "$amp_dir/TA/umi_filtered.png" \
        "$amp_dir/TA/umi_per_cell_filtered.png" || exit 1
fi
