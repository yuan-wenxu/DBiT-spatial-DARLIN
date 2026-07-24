#!/bin/bash
set -o pipefail

show_help() {
    cat << EOF
Usage: $0 <config_file>

Process image and perform cell segmentation using StarDist.

Arguments:
  config_file   Per-dataset QC configuration file

Examples:
  $0 dbit.config.sh
EOF
}

SCRIPT_DIR=${QC_SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)} || exit 1
REPO_DIR=${REPO_DIR:-$(cd "$SCRIPT_DIR/../.." && pwd)} || exit 1
PYTHON_DIR="$SCRIPT_DIR/python"

if [[ ${1:-} == -h || ${1:-} == --help ]]; then show_help; exit 0; fi
if [[ $# -ne 1 ]]; then show_help >&2; exit 1; fi
config_file=$1
if [[ ! -f "$config_file" ]]; then
    echo "Error: config file not found: $config_file" >&2; exit 1
fi

source "$config_file"
pixi_env=${image_pixi_env:-image}
pixi_env_dir=${pixi_env_dir:-$REPO_DIR}
result_path=${image_result_path:-}

if [[ -z ${image_path:-} ]]; then
    echo "Error: image_path must be set in the QC config." >&2; exit 1
fi
if [[ -z ${orientation:-} || -z ${swap_xy:-} ]]; then
    echo "Error: orientation and swap_xy must be set by the image step." >&2; exit 1
fi
case "${swap_xy,,}" in
    true) swap_xy=True ;;
    false) swap_xy=False ;;
    *) echo "Error: swap_xy must be True or False; got '$swap_xy'." >&2; exit 1 ;;
esac
for variable in x_spots_number y_spots_number length_spot interval; do
    if [[ -z ${!variable:-} ]]; then
        echo "Run this script through dbit.sh with chip already stored in the config." >&2
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

run_id=${SLURM_JOB_ID:-image_$$}
scratch_run_dir=""

cleanup_scratch() {
    local run_dir="${scratch:-}/dbit/${run_id:-}"
    if [[ -n "${run_id:-}" && -n "${scratch:-}" && -d "$run_dir" ]]; then
        rm -rf -- "$run_dir"
    fi
}

trap cleanup_scratch EXIT
trap 'exit 130' INT
trap 'exit 143' TERM HUP

run_pixi() {
  (
    cd "$pixi_env_dir" || exit 1
    pixi run -e "$pixi_env" "$@"
  )
}

# Validate inputs
if [ ! -f "$image_path" ]; then
    echo "Error: image file does not exist: $image_path" >&2
    exit 1
fi

image_dir=$(dirname "$image_path")
image_name=$(basename "$image_path")

if [ -z "$result_path" ]; then
    result_path=$image_dir
fi

result_path=$(normalize_dir_path "$result_path")
if [ -n "$scratch" ]; then
    scratch=$(normalize_dir_path "$scratch")
fi

mkdir -p "$result_path"

case "$orientation" in
    normal|horizontal|vertical|rotate)
        ;;
    *)
        echo "Error: orientation must be one of normal, horizontal, vertical, rotate" >&2
        exit 1
        ;;
esac

if [ ! -d "$pixi_env_dir" ]; then
    echo "Error: pixi environment dir does not exist: $pixi_env_dir" >&2
    exit 1
fi

if [ -n "$scratch" ]; then
    scratch_run_dir="$scratch/dbit/$run_id/image"
    mkdir -p "$scratch_run_dir"
    mkdir -p "$scratch_run_dir/result"
    cp "$image_path" "$scratch_run_dir/$image_name"
    run_image_path="$scratch_run_dir/$image_name"
    run_result_path="$scratch_run_dir/result"
else
    run_image_path="$image_path"
    run_result_path="$result_path"
fi

run_pixi python "$PYTHON_DIR/stardist_segment.py" \
  -ip "$run_image_path" \
  -r "$run_result_path" \
  --x_spots_number $x_spots_number \
  --y_spots_number $y_spots_number \
  --length_spot $length_spot \
  --interval $interval \
  --pixel_length $pixel_length \
  --put_text $put_text \
  --font_size $font_size \
  --top_value $top_value \
  --number_of_top_values $number_of_top_values \
  -m $model_name \
  -pt $prob_thresh \
  -nt $nms_thresh \
  --orientation "$orientation" \
  $([ "$swap_xy" = True ] && printf %s --swap_xy) || exit 1

run_pixi python "$PYTHON_DIR/cell_filter.py" \
  -f "$run_result_path" \
  -c $cutoff || exit 1

if [ -n "$scratch" ]; then
    cp -r "$scratch_run_dir/result"/* "$result_path/"
fi
