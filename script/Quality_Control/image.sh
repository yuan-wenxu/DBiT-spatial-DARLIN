#!/bin/bash

# Show help
show_help() {
    cat << EOF
Usage: $0 -i <image_name> -r <result_path> [OPTIONS]

Process image and perform cell segmentation using StarDist.

Required Arguments:
  -i, --image_name <path>         Path to input image file (required) (relative to result_path) (default: align.png)
  -r, --result_path <path>        Path to output result directory (required)

Image Processing Options:
  --x_spots_number <num>      Number of spots in x direction (default: 50)
  --y_spots_number <num>      Number of spots in y direction (default: 50)
  --length_spot <num>         Length of each spot in μm (default: 20)
  --interval <num>            Interval between two adjacent spots in μm (default: 20)
  --pixel_length <float>      Length of each pixel in μm (default: 0.294)
  --put_text <bool>           Whether to put text on the image (True/False, yes/no, 1/0) (default: True)
  --font_size <num>           Font size of the text (default: 1)
  --scratch <path>            Path to scratch directory for intermediate files (optional)

StarDist Quality Control Options:
  --top_value <num>              Top value threshold for image quality check (default: 50)
  --number_of_top_values <num>   Number of top values to check (default: 1500)

StarDist Detection Options:
  --model <num>                Pretrained model name (default: 2D_versatile_fluo)
  --prob_thresh <num>          Detection probability threshold (default: 0.5)
  --nms_thresh <num>           NMS IoU threshold (default: 0.6)

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
length_spot=${length_spot:-20}
interval=${interval:-20}
pixel_length=${pixel_length:-0.294}
put_text=${put_text:-True}
font_size=${font_size:-1}

# StarDist Quality Control Options
top_value=${top_value:-50}
number_of_top_values=${number_of_top_values:-1500}

# StarDist Detection Options
model_name=${model_name:-2D_versatile_fluo}
prob_thresh=${prob_thresh:-0.5}
nms_thresh=${nms_thresh:-0.6}
cutoff=${cutoff:-100}

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
        --scratch) scratch=$2; shift 2 ;;
        # StarDist Quality Control Options
        --top_value) top_value=$2; shift 2 ;;
        --number_of_top_values) number_of_top_values=$2; shift 2 ;;
        # StarDist Detection Options
        --model) model_name=$2; shift 2 ;;
        --prob_thresh) prob_thresh=$2; shift 2 ;;
        --nms_thresh) nms_thresh=$2; shift 2 ;;
        --cutoff) cutoff=$2; shift 2 ;;
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


if [ -n "$scratch" ]; then
    mkdir -p "$scratch/image"
    mkdir -p "$scratch/result"
    cp "$result_path/$image_name" "$scratch/image/$image_name"
    run_image_path="$scratch/image/$image_name"
    run_result_path="$scratch/result"
else
    run_image_path="$result_path/$image_name"
    run_result_path="$result_path"
fi

conda run -n stardist --no-capture-output python ./python/stardist_segment.py \
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
  -nt $nms_thresh

conda run -n stardist --no-capture-output python ./python/cell_filter.py \
  -f "$run_result_path" \
  -c $cutoff

if [ -n "$scratch" ]; then
    cp -r "$scratch/result"/* "$result_path/"
    rm -rf "$scratch/image"
    rm -rf "$scratch/result"
fi
