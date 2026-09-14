#!/bin/bash
set -o pipefail

SCRIPT_DIR=${STEP_SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)} || exit 1
REPO_DIR=${REPO_DIR:-$(cd "$SCRIPT_DIR/../.." && pwd)} || exit 1
PYTHON_DIR="$SCRIPT_DIR/python"
IMAGE_SCRIPT="$PYTHON_DIR/image_segment.py"
FILTER_SCRIPT="$PYTHON_DIR/image_filter.py"
MERGE_SCRIPT="$PYTHON_DIR/merge_on_image.py"

if [[ $# -ne 1 ]]; then echo "Error: expected one config file argument" >&2; exit 1; fi
config_file=$1
if [[ ! -f "$config_file" ]]; then
    echo "Error: config file not found: $config_file" >&2
    exit 1
fi

source "$config_file"
pixi_env_name=${pixi_env:-default}
pixi_env_dir=${pixi_env_dir:-$REPO_DIR}
result_path=${image_result_path:-}

if [[ -z ${image_path:-} ]]; then
    echo "Error: image_path must be set in the QC config." >&2
    exit 1
fi
for variable in x_spots_number y_spots_number length_spot interval pixel_length \
    barcode_a_whitelist_path barcode_b_whitelist_path; do
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
    local status=$?
    trap - EXIT INT TERM HUP
    if [[ -n ${scratch_run_dir:-} && -d $scratch_run_dir ]]; then
        if (( status != 0 )) && [[ -d $scratch_run_dir/result && -n ${result_path:-} ]]; then
            echo "Recovering image scratch outputs after exit status $status: $result_path" >&2
            mkdir -p "$result_path" && cp -a "$scratch_run_dir/result/." "$result_path/" || \
                echo "Warning: failed to recover scratch outputs: $scratch_run_dir/result" >&2
        fi
        rm -rf -- "$scratch_run_dir" || \
            echo "Warning: failed to clean scratch directory: $scratch_run_dir" >&2
    fi
    exit "$status"
}

enable_cleanup() {
    trap cleanup_scratch EXIT
    trap 'exit 130' INT
    trap 'exit 143' TERM
    trap 'exit 129' HUP
}

run_pixi() {
    (
        cd "$pixi_env_dir" || exit 1
        pixi run -e "$pixi_env_name" "$@"
    )
}

merge_with_image() {
    local frame
    local merge_args=(
        --image "$image_path"
        --mask "$frame_mask_path"
        --x_spots_number "$x_spots_number"
        --y_spots_number "$y_spots_number"
        --length_spot "$length_spot"
        --interval "$interval"
        --pixel_length "$pixel_length"
    )
    for frame in "$@"; do
        merge_args+=(--frame "$frame")
    done
    run_pixi python "$MERGE_SCRIPT" "${merge_args[@]}"
}

if [[ ! -f "$image_path" ]]; then
    echo "Error: image file does not exist: $image_path" >&2
    exit 1
fi
image_dir=$(dirname "$image_path")
image_name=$(basename "$image_path")
frame_mask_path="$image_dir/mask.png"
if [[ ! -f "$frame_mask_path" ]]; then
    echo "Error: adjacent frame mask does not exist: $frame_mask_path" >&2
    exit 1
fi
if [[ ! -f "$barcode_a_whitelist_path" || ! -f "$barcode_b_whitelist_path" ]]; then
    echo "Error: barcode whitelist file does not exist." >&2
    exit 1
fi

if [[ -z "$result_path" ]]; then
    result_path=$image_dir
fi
result_path=$(normalize_dir_path "$result_path")
if [[ -n ${scratch:-} ]]; then
    scratch=$(normalize_dir_path "$scratch")
fi
mkdir -p "$result_path"

if [[ ! -d "$pixi_env_dir" ]]; then
    echo "Error: pixi environment dir does not exist: $pixi_env_dir" >&2
    exit 1
fi
if [[ ! -f "$IMAGE_SCRIPT" ]]; then
    echo "Error: image processing script is missing: $IMAGE_SCRIPT" >&2
    exit 1
fi

if [[ -n ${scratch:-} ]]; then
    scratch_run_dir="$scratch/dbit/$run_id/image"
    enable_cleanup
    mkdir -p "$scratch_run_dir/result"
    cp "$image_path" "$scratch_run_dir/$image_name"
    cp "$frame_mask_path" "$scratch_run_dir/mask.png"
    run_image_path="$scratch_run_dir/$image_name"
    run_result_path="$scratch_run_dir/result"
else
    run_image_path="$image_path"
    run_result_path="$result_path"
fi

run_pixi python "$IMAGE_SCRIPT" \
    --image_path "$run_image_path" \
    --result_path "$run_result_path" \
    --barcodeA_whitelist "$barcode_a_whitelist_path" \
    --barcodeB_whitelist "$barcode_b_whitelist_path" \
    --x_spots_number "$x_spots_number" \
    --y_spots_number "$y_spots_number" \
    --length_spot "$length_spot" \
    --interval "$interval" \
    --pixel_length "$pixel_length" || exit 1

if [[ -n ${scratch:-} ]]; then
    cp -a "$scratch_run_dir/result/." "$result_path/"
fi

tissue_positions_path="$result_path/tissue_positions.tsv.gz"
grayscale_image_path="$result_path/fullres_grayscale.png"
if [[ ! -f "$grayscale_image_path" ]]; then
    echo "Error: full-resolution grayscale image is missing: $grayscale_image_path" >&2
    exit 1
fi
filter_args=(
    python "$FILTER_SCRIPT"
    --tissue_positions_file "$tissue_positions_path"
    --x_spots_number "$x_spots_number"
    --y_spots_number "$y_spots_number"
    --length_spot "$length_spot"
    --interval "$interval"
    --pixel_length "$pixel_length"
)
run_filter=false

if [[ -n ${mrna_dir:-} ]]; then
    mrna_dir=$(normalize_dir_path "$mrna_dir")
    filter_args+=(--mrna_path "$mrna_dir")
    run_filter=true
fi

if [[ -n ${darlin_dir:-} ]]; then
    darlin_dir=$(normalize_dir_path "$darlin_dir")
    filter_args+=(--darlin_path "$darlin_dir")
    run_filter=true
fi

if $run_filter; then
    if [[ -n ${mrna_dir:-} ]]; then
        run_pixi "${filter_args[@]}" \
            2>&1 | tee "$mrna_dir/filtered_plot.log" || exit 1
    else
        run_pixi "${filter_args[@]}" || exit 1
    fi
fi

if [[ -n ${mrna_dir:-} ]]; then
    if [[ -n ${mrna_output_path:-} ]]; then
        mrna_matrix_dir="$(normalize_dir_path "$mrna_output_path")/matrix"
    else
        mrna_matrix_dir=$(realpath -m "$mrna_dir/../../../matrix")
    fi
    mkdir -p "$mrna_matrix_dir" || exit 1
    cp -f -- "$tissue_positions_path" "$mrna_matrix_dir/tissue_positions.tsv.gz" || exit 1
    cp -f -- "$grayscale_image_path" "$mrna_matrix_dir/tissue_raw_image.png" || exit 1
    echo "Copied tissue positions: $mrna_matrix_dir/tissue_positions.tsv.gz"
    echo "Copied grayscale image: $mrna_matrix_dir/tissue_raw_image.png"
fi

if [[ -n ${mrna_dir:-} ]]; then
    merge_with_image \
        "$mrna_dir/raw/umap_filtered.png" \
        "$mrna_dir/raw/umi_filtered.png" \
        "$mrna_dir/raw/gene_filtered.png" || exit 1
fi
