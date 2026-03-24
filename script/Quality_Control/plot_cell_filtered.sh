#!/bin/bash

# Show help message
show_help() {
    cat << EOF
Usage: $0 -c <cell_number_file> [ -d <mrna_dir> ] [ -a <amp_dir> -w <whitelist> ] [OPTIONS]

Plot filtered results (cell-filtered) for mRNA and/or amplicon data.

Required Arguments:
  -c, --cell_number_file <path>     Path to cell number file (e.g. cell_num_area.csv)
  -w, --whitelist <path>            Path to barcode whitelist file
  
Optional (at least one of the two groups below):
  -d, --mrna_dir <path>             Path to mRNA results directory
  -a, --amp_dir <path>              Path to amplicon results directory

Spatial/Plot Options (passed to both mrna_cell.py and amp_cell.py):
  --x_spots_number <num>            Number of spots in x direction (default: 50)
  --y_spots_number <num>            Number of spots in y direction (default: 50)
  --length_spot <num>               Length of each spot in pixels (default: 50)
  --interval <num>                  Interval between spots in pixels (default: 50)
  --pixel_length <float>            Length of each pixel in microns (default: 0.294)

Other Options:
  -h, --help                        Show this help message and exit

Examples:
  # Both mRNA and amplicon
  $0 -c cell_num_area.csv -d /path/to/mrna/results -a /path/to/amp/results -w barcodes.tsv

EOF
}

# Set default values (order matches the help message)

MERGE_SCRIPT=./python/image_process/merge_on_gray.py

merge_with_gray() {
    local search_dir="$1"
    local gray="$(dirname "$cell_number_file")/gray.png"
    if [[ -f "$gray" ]]; then
        python "$MERGE_SCRIPT" \
            --gray "$gray" \
            --search-dir "$search_dir" \
            --recursive
    else
        echo "Warning: gray.png not found at $gray, skipping merge."
    fi
}

# Spatial/Plot Options
x_spots_number=${x_spots_number:-50}
y_spots_number=${y_spots_number:-50}
length_spot=${length_spot:-50}
interval=${interval:-50}
pixel_length=${pixel_length:-0.294}

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
while getopts "c:d:a:b:w:h" opt; do
    case $opt in
        c) cell_number_file=$OPTARG ;;
        d) mrna_dir=$OPTARG ;;
        a) amp_dir=$OPTARG ;;
        b) darlin_from_bam_dir=$OPTARG ;;
        w) whitelist_path=$OPTARG ;;
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
        --x_spots_number) x_spots_number=$2; shift 2 ;;
        --y_spots_number) y_spots_number=$2; shift 2 ;;
        --length_spot) length_spot=$2; shift 2 ;;
        --interval) interval=$2; shift 2 ;;
        --pixel_length) pixel_length=$2; shift 2 ;;
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
if [ -z "$whitelist_path" ]; then
    echo "Error: -w (whitelist) is required" >&2
    exit 1
fi
if [ -z "$mrna_dir" ] && [ -z "$amp_dir" ]; then
    echo "Error: at least one of -d (mrna_dir) or -a (amp_dir) must be specified" >&2
    echo "Use -h or --help for usage information" >&2
    exit 1
fi

# mRNA cell-filtered plot (only when -d is provided)
if [ -n "$mrna_dir" ]; then
    python ./python/mrna_cell.py \
        -c "$cell_number_file" \
        -d "$mrna_dir" \
        --x_spots_number "$x_spots_number" \
        --y_spots_number "$y_spots_number" \
        --length_spot "$length_spot" \
        --interval "$interval" \
        --pixel_length "$pixel_length" \
        &> "$mrna_dir/filtered_plot.log"
    merge_with_gray "$mrna_dir"
fi

# Amplicon cell-filtered plot (only when -a is provided)
if [ -n "$amp_dir" ]; then
    python ./python/amplicon_cell.py \
        -c "$cell_number_file" \
        -d "$amp_dir" \
        -w "$whitelist_path" \
        --x_spots_number "$x_spots_number" \
        --y_spots_number "$y_spots_number" \
        --length_spot "$length_spot" \
        --interval "$interval" \
        --pixel_length "$pixel_length" \
        &> "$amp_dir/filtered_plot.log"
    merge_with_gray "$amp_dir"
fi
