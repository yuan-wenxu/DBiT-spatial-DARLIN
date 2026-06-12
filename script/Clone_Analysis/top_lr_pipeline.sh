#!/bin/bash

set -euo pipefail

show_help() {
    cat << EOF
Usage: $0 -i <input_dir> -b <bank_dir> --cluster-csv <csv> [options]

Run the top-LR plotting pipeline in this order:
  1. cellcount_filter.py
  2. allele_bank_filter.py
  3. top_lr_plot.py

Required Arguments:
  -i, --input-dir <dir>             Path containing CA/RA/TA subdirectories
  -b, --bank-dir <dir>              Allele bank directory
      --cluster-csv <csv>           mRNA cluster CSV used as the plotting background

Other Options:
      --labels <label...>           Labels to process (default: CA RA TA)
      --output-dir <dir>            Output directory (default: input_dir)
      --min-sequence-length <int>   Minimum sequence length for allele analysis
      --top-n <int>                 Number of top LR plots per label (default: 10)
      --rotate <0|90|180|270>       Counterclockwise plot rotation (default: 0)
  -h, --help                        Show this help message and exit

EOF
}


labels=("CA" "RA" "TA")
output_dir=""
min_sequence_length=20
top_n=10
rotate=0

cellcount_filter="./python/cellcount_filter.py"
allele_bank_filter="./python/allele_bank_filter.py"
top_lr_plot="./python/top_lr_plot.py"


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
while getopts "i:b:h" opt; do
    case $opt in
        i) input_dir=$OPTARG ;;
        b) bank_dir=$OPTARG ;;
        h) show_help; exit 0 ;;
        ?) echo "Invalid option: -$OPTARG" >&2
           echo "Use -h or --help for usage information" >&2
           exit 1 ;;
    esac
done


set -- "${long_args[@]}"
while [[ $# -gt 0 ]]; do
    case $1 in
        --input-dir) input_dir=$2; shift 2 ;;
        --bank-dir) bank_dir=$2; shift 2 ;;
        --cluster-csv) cluster_csv=$2; shift 2 ;;
        --output-dir) output_dir=$2; shift 2 ;;
        --min-sequence-length) min_sequence_length=$2; shift 2 ;;
        --top-n) top_n=$2; shift 2 ;;
        --rotate) rotate=$2; shift 2 ;;
        --labels)
            labels=()
            shift
            while [[ $# -gt 0 ]] && [[ "$1" != --* ]]; do
                labels+=("$1")
                shift
            done
            ;;
        --help) show_help; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done


if [ -z "${input_dir:-}" ]; then
    echo "Error: -i (--input-dir) is required" >&2
    exit 1
fi

if [ -z "${bank_dir:-}" ]; then
    echo "Error: -b (--bank-dir) is required" >&2
    exit 1
fi

if [ -z "${cluster_csv:-}" ]; then
    echo "Error: --cluster-csv is required" >&2
    exit 1
fi


run_step() {
    local title=$1
    shift

    echo
    echo "=== $title ==="
    python "$@"
}


cellcount_args=(--input-dir "$input_dir" --cluster-csv "$cluster_csv" --labels "${labels[@]}")
allele_args=(--input_dir "$input_dir" --bank-dir "$bank_dir" --min-sequence-length "$min_sequence_length" --labels "${labels[@]}")
plot_args=(--input-dir "$input_dir" --cluster-csv "$cluster_csv" --top-n "$top_n" --rotate "$rotate" --labels "${labels[@]}")

if [ -n "$output_dir" ]; then
    cellcount_args+=(--output-dir "$output_dir")
    allele_args+=(--output-dir "$output_dir")
    plot_args+=(--output-dir "$output_dir")
fi


run_step \
    "Step 1/3: Cell Count Split" \
    "$cellcount_filter" \
    "${cellcount_args[@]}"

run_step \
    "Step 2/3: Allele Bank Filter" \
    "$allele_bank_filter" \
    "${allele_args[@]}"

run_step \
    "Step 3/3: Top LR Plot" \
    "$top_lr_plot" \
    "${plot_args[@]}"

echo
echo "Pipeline completed successfully."
