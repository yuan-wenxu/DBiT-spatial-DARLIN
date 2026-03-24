#!/bin/bash

# Show help
show_help() {
    cat << EOF
Usage: $0 -i <image_path> -r <result_path> [OPTIONS]

Process image and perform cell segmentation using Cellpose.

Required Arguments:
  -i, --image_name <path>         Path to input image file (required) (relative to result_path) (default: align.png)
  -r, --result_path <path>        Path to output result directory (required)

Image Processing Options:
  --x_spots_number <num>      Number of spots in x direction (default: 50)
  --y_spots_number <num>      Number of spots in y direction (default: 50)
  --length_spot <num>         Length of each spot in μm (default: 50)
  --interval <num>            Interval between two adjacent spots in μm (default: 50)
  --pixel_length <float>      Length of each pixel in μm (default: 0.294)
  --put_text <bool>           Whether to put text on the image (True/False, yes/no, 1/0) (default: True)
  --font_size <num>           Font size of the text (default: 1)

Cellpose Quality Control Options:
  --top_value <num>              Top value threshold for image quality check (default: 50)
  --number_of_top_values <num>   Number of top values to check (default: 1500)

Cellpose Detection Options:
  --dmin <num>                Minimum cell diameter in pixels (default: 20)
  --dmax <num>                Maximum cell diameter in pixels (default: 50)
  --step <num>                Step size for diameter search (default: 10)
  --cutoff <num>              Minimum size for cell detection (default: 75)

Cellpose Image Processing Options:
  --photo_size <num>          Image tile size for processing in pixels (default: 170)
  --photo_step <num>          Step size for image tiling in pixels (default: 170)

Other Options:
  -h, --help                      Show this help message and exit

Examples:
  # Basic usage with default values
  $0 -i image.tif -r /path/to/result

EOF
}

# Set default values
image_name=${image_name:-align.png}
# Image Processing Options
x_spots_number=${x_spots_number:-50}
y_spots_number=${y_spots_number:-50}
length_spot=${length_spot:-50}
interval=${interval:-50}
pixel_length=${pixel_length:-0.294}
put_text=${put_text:-True}
font_size=${font_size:-1}

# Cellpose Quality Control Options
top_value=${top_value:-50}
number_of_top_values=${number_of_top_values:-1500}

# Cellpose Detection Options
dmin=${dmin:-20}
dmax=${dmax:-50}
step=${step:-10}
cutoff=${cutoff:-75}

# Cellpose Image Processing Options
photo_size=${photo_size:-170}
photo_step=${photo_step:-170}

# Short options
short_args=()
while [[ $# -gt 0 ]]; do
    case $1 in
        --*) break ;;
        *) short_args+=("$1"); shift ;;
    esac
done
long_args=("$@")

set -- "${short_args[@]}"
OPTIND=1

while getopts "i:r:h" opt; do
    case $opt in
        i) image_name=$OPTARG ;;
        r) result_path=$OPTARG ;;
        h) show_help; exit 0 ;;
        ?) echo "Invalid option: -$OPTARG" >&2
            echo "Use -h or --help for usage information"
            exit 1 ;;
    esac
done

# Long options
set -- "${long_args[@]}"
while [[ $# -gt 0 ]]; do
    case $1 in
        # Image Processing Options
        --x_spots_number) x_spots_number=$2; shift 2 ;;
        --y_spots_number) y_spots_number=$2; shift 2 ;;
        --pixel_length) pixel_length=$2; shift 2 ;;
        --length_spot) length_spot=$2; shift 2 ;;
        --interval) interval=$2; shift 2 ;;
        --pixel_length) pixel_length=$2; shift 2 ;;
        --put_text) put_text=$2; shift 2 ;;
        --font_size) font_size=$2; shift 2 ;;
        # Cellpose Quality Control Options
        --top_value) top_value=$2; shift 2 ;;
        --number_of_top_values) number_of_top_values=$2; shift 2 ;;
        # Cellpose Detection Options
        --dmin) dmin=$2; shift 2 ;;
        --dmax) dmax=$2; shift 2 ;;
        --step) step=$2; shift 2 ;;
        # Cellpose Image Processing Options
        --photo_size) photo_size=$2; shift 2 ;;
        --photo_step) photo_step=$2; shift 2 ;;
        # Other Options
        --help) show_help; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

# Check required arguments
if [ -z "$image_name" ] || [ -z "$result_path" ]; then
    echo "Error: -i (image_name) and -r (result_path) are required" >&2
    echo "Use -h or --help for usage information" >&2
    exit 1
fi


conda run -n py38 python ./python/image.py \
    -ip "$result_path/$image_name" \
    -r "$result_path" \
    -x "$x_spots_number" \
    -y "$y_spots_number" \
    -p "$pixel_length" \
    -l "$length_spot" \
    -i "$interval" \
    -t "$put_text" \
    -fs "$font_size" \
    -top_value "$top_value" \
    -number_of_top_values "$number_of_top_values" \
    -dmin "$dmin" \
    -dmax "$dmax" \
    -step "$step" \
    -photo_size "$photo_size" \
    -photo_step "$photo_step"

conda run -n py38 python ./python/cell_filter.py \
    -f "$result_path" \
    -c "$cutoff"
