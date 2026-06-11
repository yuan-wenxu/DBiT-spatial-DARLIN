#!/bin/bash

# Show help message
show_help() {
    cat << EOF
Usage: $0 -c <cell_number_file> [ -m <mrna_dir> ] [ -a <amp_dir> -w <whitelist> ] [OPTIONS]

Plot filtered results (cell-filtered) for mRNA and/or amplicon data.

Required Arguments:
  -c, --cell_number_file <path>     Path to cell number file (e.g. filtered_results.csv)
  
Optional (at least one of the two groups below):
  -m, --mrna_dir <path>             Path to mRNA results directory
  -a, --amp_dir <path>              Path to amplicon results directory
  -w, --whitelist <path>            Path to barcode whitelist file (required only when -a/--amp_dir is provided)
  -g, --gray_path <path>            Path to gray image for merge (optional; default: gray.png next to cell_number_file)
  --orientation <mode>          Transform filtered images before merging: normal, horizontal, vertical, rotate (default: normal)
  --swap_xy                         Swap x and y axes after applying orientation before merging (optional)

Spatial/Plot Options (passed to both mrna_cell.py and amp_cell.py):
  --x_spots_number <num>            Number of spots in x direction (default: 50)
  --y_spots_number <num>            Number of spots in y direction (default: 50)
  --length_spot <num>               Length of each spot in pixels (default: 20)
  --interval <num>                  Interval between spots in pixels (default: 20)
  --pixel_length <float>            Length of each pixel in microns (default: 0.294)

Orientation Notes:
  normal                            Keep the coordinate system unchanged
  horizontal                        Flip the coordinate system horizontally
  vertical                          Flip the coordinate system vertically
  rotate                            Rotate the coordinate system 180 degrees; equivalent to horizontal + vertical
  --swap_xy                         Swap x and y coordinate axes after applying orientation
  horizontal + --swap_xy            Rotate the coordinate system 90 degrees counterclockwise
  vertical + --swap_xy              Rotate the coordinate system 90 degrees clockwise

Pixi environment options:
  --pixi_env <name>                   Name of the Pixi environment to use (optional; default: dbit)
  --pixi_env_dir <path>               Directory containing pixi.toml (optional; default: repository root)

Other Options:
  -h, --help                        Show this help message and exit

Examples:
  # mRNA only
  $0 -c filtered_results.csv -m /path/to/mrna/results

  # Use a custom gray image for merge
  $0 -c filtered_results.csv -m /path/to/mrna/results -g /path/to/gray.png

  # Rotate filtered images before merging on gray
  $0 -c filtered_results.csv -m /path/to/mrna/results -o rotate

  # Both mRNA and amplicon
  $0 -c filtered_results.csv -m /path/to/mrna/results -a /path/to/amp/results -w barcodes.tsv

EOF
}

# Set default values (order matches the help message)

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd) || exit 1
PYTHON_DIR="$SCRIPT_DIR/python"
MERGE_SCRIPT="$PYTHON_DIR/image_process/merge_on_gray.py"
pixi_env_dir=${pixi_env_dir:-$(cd "$SCRIPT_DIR/../.." && pwd)}

normalize_dir_path() {
    local path="$1"
    while [[ "$path" != "/" && "$path" == */ ]]; do
        path="${path%/}"
    done
    printf '%s\n' "$path"
}

merge_with_gray() {
    local search_dir="$1"
    if [[ -f "$gray_path" ]]; then
        (
            cd "$pixi_env_dir" || exit 1
            merge_args=(
                --gray "$gray_path"
                --search-dir "$search_dir"
                --orientation "$orientation"
                --recursive
            )
            if [[ "$swap_xy" == True ]]; then
                merge_args+=(--swap_xy)
            fi
            pixi run -e "$pixi_env" python "$MERGE_SCRIPT" "${merge_args[@]}"
        )
    else
        echo "Warning: gray image not found at $gray_path, skipping merge."
    fi
}

# Spatial/Plot Options
x_spots_number=${x_spots_number:-50}
y_spots_number=${y_spots_number:-50}
length_spot=${length_spot:-20}
interval=${interval:-20}
pixel_length=${pixel_length:-0.294}
orientation=${orientation:-normal}
swap_xy=${swap_xy:-False}
pixi_env=${pixi_env:-dbit}

# Parse only short options with getopts; stop on --/long options to avoid getopts treating "--" as invalid.
short_args=()
while [[ $# -gt 0 ]]; do
    case $1 in
        --*) break ;;
        *) short_args+=("$1"); shift ;;
    esac
done
long_args=("$@")

# Parse short options
set -- "${short_args[@]}"
OPTIND=1
while getopts "c:m:a:w:g:h" opt; do
    case $opt in
        c) cell_number_file=$OPTARG ;;
        m) mrna_dir=$OPTARG ;;
        a) amp_dir=$OPTARG ;;
        w) whitelist_path=$OPTARG ;;
        g) gray_path=$OPTARG ;;
        h) show_help; exit 0 ;;
        ?) echo "Invalid option: -$OPTARG" >&2
            echo "Use -h or --help for usage information" >&2
            exit 1 ;;
    esac
done

# Parse long options
set -- "${long_args[@]}"
while [[ $# -gt 0 ]]; do
    case $1 in
        --cell_number_file) cell_number_file=$2; shift 2 ;;
        --mrna_dir) mrna_dir=$2; shift 2 ;;
        --amp_dir) amp_dir=$2; shift 2 ;;
        --whitelist) whitelist_path=$2; shift 2 ;;
        --gray_path) gray_path=$2; shift 2 ;;
        --orientation) orientation=$2; shift 2 ;;
        --swap_xy) swap_xy=True; shift ;;
        --x_spots_number) x_spots_number=$2; shift 2 ;;
        --y_spots_number) y_spots_number=$2; shift 2 ;;
        --length_spot) length_spot=$2; shift 2 ;;
        --interval) interval=$2; shift 2 ;;
        --pixel_length) pixel_length=$2; shift 2 ;;
        --pixi_env) pixi_env=$2; shift 2 ;;
        --pixi_env_dir) pixi_env_dir=$2; shift 2 ;;
        --help) show_help; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

# Validate required arguments: cell_number_file is required, and at least one data directory must be provided.
if [ -z "$cell_number_file" ]; then
    echo "Error: -c (cell_number_file) is required" >&2
    echo "Use -h or --help for usage information" >&2
    exit 1
fi
if [ -z "$mrna_dir" ] && [ -z "$amp_dir" ]; then
    echo "Error: at least one of -m (mrna_dir) or -a (amp_dir) must be specified" >&2
    echo "Use -h or --help for usage information" >&2
    exit 1
fi
if [ -n "$amp_dir" ] && [ -z "$whitelist_path" ]; then
    echo "Error: -w (whitelist) is required when -a (amp_dir) is specified" >&2
    echo "Use -h or --help for usage information" >&2
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
        echo "Error: -o/--orientation must be one of normal, horizontal, vertical, rotate" >&2
        exit 1
        ;;
esac

if [[ ! -d "$pixi_env_dir" ]]; then
    echo "Error: pixi environment dir does not exist: $pixi_env_dir" >&2
    exit 1
fi

# mRNA cell-filtered plot (only when -m is provided)
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
    ) &> "$mrna_dir/filtered_plot.log"
    merge_with_gray "$mrna_dir"
fi

# Amplicon cell-filtered plot (only when -a is provided)
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
    ) &> "$amp_dir/filtered_plot.log"
    merge_with_gray "$amp_dir"
fi
