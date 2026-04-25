#!/bin/bash

set -euo pipefail

# Show help message
show_help() {
	cat << EOF
Usage: $0 -i <input_dir> --labels RA TA

Run the clone analysis pipeline in this order:
  1. allele_bank_filter.py
  2. cellcount_filter.py
  3. distance_filter.py
  4. clustered_clone_plot.py

Required Arguments:
  -i, --input_dir <dir>             Path containing CA/RA/TA subdirectories

Other Options:
	--labels <label...>               Labels to process (default: CA RA TA)
	--output-dir <dir>                Output directory (default: input_dir)
	--bank-dir <dir>                  Allele bank directory
	--min-sequence-length <int>       Minimum sequence length for allele analysis
	--force-reanalyze                 Reanalyze sequences even if cache exists
	--distance-unit <spot|um>         Distance unit for distance analysis
	--length-spot <num>               Spot size in um
	--interval <num>                  Center-to-center gap contribution in um
	--cluster-threshold <num>         Cluster threshold for mean pairwise distance
	--summary-subdir <name>           Distance summary subdirectory name
	--pixel-length <num>              Pixel size for clustered clone plotting
	--x-spots-number <int>            Number of spots in x direction
	--y-spots-number <int>            Number of spots in y direction
  -h, --help                          Show this help message and exit

EOF
}


# Set default values
labels=("CA" "RA" "TA")
output_dir=""
bank_dir="/mnt/dbit/data/reference/allele_bank"
min_sequence_length=20
force_reanalyze=false
distance_unit="um"
length_spot=20
interval=20
cluster_threshold=200
summary_subdir="distance_analysis"
pixel_length=0.294
x_spots_number=50
y_spots_number=50

allele_bank_filter="./python/allele_bank_filter.py"
cellcount_filter="./python/cellcount_filter.py"
distance_filter="./python/distance_filter.py"
clustered_clone_plot="./python/clustered_clone_plot.py"


# Parse only short options first; stop when long options begin.
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
while getopts "i:h" opt; do
	case $opt in
		i) input_dir=$OPTARG ;;
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
		--input_dir) input_dir=$2; shift 2 ;;
		--output-dir) output_dir=$2; shift 2 ;;
		--bank-dir) bank_dir=$2; shift 2 ;;
		--min-sequence-length) min_sequence_length=$2; shift 2 ;;
		--force-reanalyze) force_reanalyze=true; shift ;;
		--labels)
			labels=()
			shift
			while [[ $# -gt 0 ]] && [[ "$1" != --* ]]; do
				labels+=("$1")
				shift
			done
			;;
		--distance-unit) distance_unit=$2; shift 2 ;;
		--length-spot) length_spot=$2; shift 2 ;;
		--interval) interval=$2; shift 2 ;;
		--cluster-threshold) cluster_threshold=$2; shift 2 ;;
		--summary-subdir) summary_subdir=$2; shift 2 ;;
		--pixel-length) pixel_length=$2; shift 2 ;;
		--x-spots-number) x_spots_number=$2; shift 2 ;;
		--y-spots-number) y_spots_number=$2; shift 2 ;;
		--help) show_help; exit 0 ;;
		*) echo "Unknown option: $1" >&2; exit 1 ;;
	esac
done


# Check required arguments
if [ -z "$input_dir" ]; then
	echo "Error: -i (--input_dir) is required" >&2
	echo "Use -h or --help for usage information" >&2
	exit 1
fi

run_step() {
	local title=$1
	shift

	echo
	echo "=== $title ==="
	python "$@"
}


common_args=(--input_dir "$input_dir" --labels "${labels[@]}")
allele_args=("${common_args[@]}" --bank-dir "$bank_dir" --min-sequence-length "$min_sequence_length")
cellcount_args=("${common_args[@]}")
distance_args=("${common_args[@]}" --distance-unit "$distance_unit" --length-spot "$length_spot" --interval "$interval" --cluster-threshold "$cluster_threshold")
clustered_plot_args=("${common_args[@]}" --summary-subdir "$summary_subdir" --length-spot "$length_spot" --interval "$interval" --pixel-length "$pixel_length" --x-spots-number "$x_spots_number" --y-spots-number "$y_spots_number")

if [ -n "$output_dir" ]; then
	allele_args+=(--output-dir "$output_dir")
	cellcount_args+=(--output-dir "$output_dir")
	distance_args+=(--output-dir "$output_dir")
	clustered_plot_args+=(--output-dir "$output_dir")
fi

if [ "$force_reanalyze" = true ]; then
	allele_args+=(--force-reanalyze)
fi


run_step \
	"Step 1/4: Allele Bank Filter" \
	"$allele_bank_filter" \
	"${allele_args[@]}"

run_step \
	"Step 2/4: Cell Count Filter" \
	"$cellcount_filter" \
	"${cellcount_args[@]}"

run_step \
	"Step 3/4: Distance Filter" \
	"$distance_filter" \
	"${distance_args[@]}"

run_step \
	"Step 4/4: Clustered Clone Plot" \
	"$clustered_clone_plot" \
	"${clustered_plot_args[@]}"

echo
echo "Pipeline completed successfully."
