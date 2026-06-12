#!/bin/bash

# Show help
show_help() {
    cat << EOF
Usage: $0 -i <image_path> [-r <result_path>] [-o <orientation>] [OPTIONS]

Process image and perform cell segmentation using StarDist.

Required Arguments:
  -i, --image_path <path>         Path to input image file (required)

Output Options:
  -r, --result_path <path>        Path to output result directory (optional; default: image file directory)

Image Processing Options:
  --x_spots_number <num>      Number of spots in x direction (default: 50)
  --y_spots_number <num>      Number of spots in y direction (default: 50)
  --length_spot <num>         Length of each spot in μm (default: 20)
  --interval <num>            Interval between two adjacent spots in μm (default: 20)
  --pixel_length <float>      Length of each pixel in μm (default: 0.294)
  --put_text <bool>           Whether to put text on the image (True/False, yes/no, 1/0) (default: True)
  --font_size <num>           Font size of the text (default: 1)
  --orientation <mode>    Grid origin orientation: normal, horizontal, vertical, rotate (default: normal)
  --swap_xy                   Swap x and y grid axes after applying orientation (optional)
  --scratch <path>            Path to scratch directory for intermediate files (optional)

Orientation Notes:
  normal                      Keep the coordinate system unchanged
  horizontal                  Flip the coordinate system horizontally
  vertical                    Flip the coordinate system vertically
  rotate                      Rotate the coordinate system 180 degrees; equivalent to horizontal + vertical
  --swap_xy                   Swap x and y coordinate axes
  horizontal + --swap_xy      Rotate the coordinate system 90 degrees counterclockwise
  vertical + --swap_xy        Rotate the coordinate system 90 degrees clockwise

StarDist Quality Control Options:
  --top_value <num>              Top value threshold for image quality check (default: 50)
  --number_of_top_values <num>   Number of top values to check (default: 1500)

StarDist Detection Options:
  --model <num>                Pretrained model name (default: 2D_versatile_fluo)
  --prob_thresh <num>          Detection probability threshold (default: 0.5)
  --nms_thresh <num>           NMS IoU threshold (default: 0.6)

Pixi environment options:
  --pixi_env <name>                   Name of the Pixi environment to use (optional; default: stardist)
  --pixi_env_dir <path>               Directory containing pixi.toml (optional; default: repository root)

Other Options:
  -h, --help                      Show this help message and exit

Examples:
  # Basic usage with default values
  $0 -i /path/to/image.tif

  # Write results to a different directory
  $0 -i /path/to/image.tif -r /path/to/result

  # Make 0_0 start at the right-bottom corner
  $0 -i /path/to/image.tif -o rotate

EOF
}

# Set default values
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd) || exit 1
PYTHON_DIR="$SCRIPT_DIR/python"

normalize_dir_path() {
    local path="$1"
    while [[ "$path" != "/" && "$path" == */ ]]; do
        path="${path%/}"
    done
    printf '%s\n' "$path"
}

# Image Processing Options
x_spots_number=${x_spots_number:-50}
y_spots_number=${y_spots_number:-50}
length_spot=${length_spot:-20}
interval=${interval:-20}
pixel_length=${pixel_length:-0.294}
put_text=${put_text:-True}
font_size=${font_size:-1}
orientation=${orientation:-normal}
swap_xy=${swap_xy:-False}

# StarDist Quality Control Options
top_value=${top_value:-50}
number_of_top_values=${number_of_top_values:-1500}

# StarDist Detection Options
model_name=${model_name:-2D_versatile_fluo}
prob_thresh=${prob_thresh:-0.5}
nms_thresh=${nms_thresh:-0.6}
cutoff=${cutoff:-100}

# Pixi environment options
pixi_env=${pixi_env:-stardist}
pixi_env_dir=${pixi_env_dir:-$(cd "$SCRIPT_DIR/../.." && pwd)}

run_pixi() {
  (
    cd "$pixi_env_dir" || exit 1
    pixi run -e "$pixi_env" "$@"
  )
}

while [[ $# -gt 0 ]]; do
    case $1 in
        # Required Arguments
        -i|--image_path)
            image_path=$2
            shift 2
            ;;
        # Output Options
        -r|--result_path)
            result_path=$2
            shift 2
            ;;
        # Image Processing Options
        --x_spots_number) x_spots_number=$2; shift 2 ;;
        --y_spots_number) y_spots_number=$2; shift 2 ;;
        --pixel_length) pixel_length=$2; shift 2 ;;
        --length_spot) length_spot=$2; shift 2 ;;
        --interval) interval=$2; shift 2 ;;
        --put_text) put_text=$2; shift 2 ;;
        --font_size) font_size=$2; shift 2 ;;
        --orientation) orientation=$2; shift 2 ;;
        --swap_xy) swap_xy=True; shift ;;
        --scratch) scratch=$2; shift 2 ;;
        # StarDist Quality Control Options
        --top_value) top_value=$2; shift 2 ;;
        --number_of_top_values) number_of_top_values=$2; shift 2 ;;
        # StarDist Detection Options
        --model) model_name=$2; shift 2 ;;
        --prob_thresh) prob_thresh=$2; shift 2 ;;
        --nms_thresh) nms_thresh=$2; shift 2 ;;
        --cutoff) cutoff=$2; shift 2 ;;
        --pixi_env) pixi_env=$2; shift 2 ;;
        --pixi_env_dir) pixi_env_dir=$2; shift 2 ;;
        # Other Options
        -h|--help) show_help; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

# Check required arguments
if [ -z "$image_path" ]; then
    echo "Error: -i/--image_path is required" >&2
    echo "Use -h or --help for usage information" >&2
    exit 1
fi

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
        echo "Error: --orientation must be one of normal, horizontal, vertical, rotate" >&2
        exit 1
        ;;
esac

if [ ! -d "$pixi_env_dir" ]; then
    echo "Error: pixi environment dir does not exist: $pixi_env_dir" >&2
    exit 1
fi

if [ -n "$scratch" ]; then
    mkdir -p "$scratch/image"
    mkdir -p "$scratch/image/result"
    cp "$image_path" "$scratch/image/$image_name"
    run_image_path="$scratch/image/$image_name"
    run_result_path="$scratch/image/result"
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
  $([ "$swap_xy" = True ] && printf %s --swap_xy)

run_pixi python "$PYTHON_DIR/cell_filter.py" \
  -f "$run_result_path" \
  -c $cutoff

if [ -n "$scratch" ]; then
    cp -r "$scratch/image/result"/* "$result_path/"
    rm -rf "$scratch/image"
fi
