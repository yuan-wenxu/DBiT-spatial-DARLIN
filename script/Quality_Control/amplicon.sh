#!/bin/bash

show_help() {
    cat << EOF
Usage: $0 <config_file>

Process amplification sequencing data with preprocessing, DARLIN correction, and visualization.

Arguments:
  config_file   Per-dataset QC configuration file

Examples:
  $0 config.sh
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
pixi_env=${pixi_env:-default}
pixi_env_dir=${pixi_env_dir:-$REPO_DIR}
initial_reads_cutoff=${initial_reads_cutoff:-100}
major_fraction_threshold_molecule=${major_fraction_threshold_molecule:-0.8}
reads_fraction_mode=${reads_fraction_mode:-sum}
reads_cutoff=${reads_cutoff:-10}
slope_cutoff=${slope_cutoff:-10}
fastq_path=$amplicon_fastq_path
output_path=${amplicon_output_path:-}
cores=${amp_cores}
cutadapt=${cutadapt}

if [[ -z ${amplicon_fastq_path:-} ]]; then
    echo "Error: amplicon_fastq_path must be set in the QC config." >&2; exit 1
fi
if [[ -z ${whitelist_path:-} ]]; then
    echo "Run this script through dbit.sh so --chip is resolved." >&2
    exit 1
fi

normalize_dir_path() {
    local path="$1"
    while [[ "$path" != "/" && "$path" == */ ]]; do
        path="${path%/}"
    done
    printf '%s\n' "$path"
}

run_pixi() {
    (
        cd "$pixi_env_dir" || exit 1
        pixi run -e "$pixi_env" "$@"
    )
}

compress_fastq_file() {
    local fq="$1"
    local threads="$2"
    local level="$3"
    if command -v pigz >/dev/null 2>&1; then
        pigz -f -p "$threads" "-$level" "$fq"
    else
        gzip -f "-$level" "$fq"
    fi
}

# Validate inputs
if [ -z "$whitelist_path" ]; then
    echo "Error: whitelist_path is required in config" >&2
    exit 1
fi

if [[ "$reads_fraction_mode" != "sum" && "$reads_fraction_mode" != "max" ]]; then
    echo "Error: reads_fraction_mode must be 'sum' or 'max'" >&2
    exit 1
fi

fastq_path=$(normalize_dir_path "$fastq_path")
if [ -n "$output_path" ]; then
    output_path=$(normalize_dir_path "$output_path")
fi
if [ -n "$scratch" ]; then
    scratch=$(normalize_dir_path "$scratch")
fi

if [ ! -d "$pixi_env_dir" ]; then
    echo "Error: pixi environment dir does not exist: $pixi_env_dir" >&2
    exit 1
fi

if [ ! -d "$fastq_path" ]; then
    echo "Error: fastq directory does not exist: $fastq_path" >&2
    exit 1
fi

if [ -z "$output_path" ]; then
    output_path=$(dirname "$fastq_path")
fi

mkdir -p "$output_path"

if [ -n "$scratch" ]; then
    scratch_input="$scratch/amplicon/input"
    scratch_output="$scratch/amplicon/output"
    mkdir -p "$scratch_input" "$scratch_output"
    cp -r "$fastq_path"/* "$scratch_input/"
    orig_output_path="$output_path"
    file_path="$scratch_input"
    output_path="$scratch_output"
else
    file_path="$fastq_path"
fi

for r1 in "$file_path"/*_R1.fq.gz; do
    [ -e "$r1" ] || { echo "Error: no *_R1.fq.gz files found in $file_path" >&2; exit 1; }
    sample_name=$(basename $r1 | sed 's/_R1.fq.gz//')
    locus=$(echo "$sample_name" | grep -oE 'CA|RA|TA')
    r2=$file_path/$sample_name"_R2.fq.gz"
    nonlocus_sample_name=$(echo "$sample_name" | sed 's/\(-CA\|-RA\|-TA\|_CA\|_RA\|_TA\)//')

    # Cutadapt and extract UMI and barcode
    gzip_after_enabled=false
    if [[ "${gzip_after_preprocess,,}" =~ ^(true|yes|1)$ ]]; then
        gzip_after_enabled=true
    fi

    if [[ "${gzip_output,,}" =~ ^(false|no|0)$ ]]; then
        preprocess_bc_ext="fq"
    else
        preprocess_bc_ext="fq.gz"
    fi
    bc_ext="$preprocess_bc_ext"
    if [ "$preprocess_bc_ext" = "fq" ] && $gzip_after_enabled; then
        bc_ext="fq.gz"
    fi

    run_pixi python "$PYTHON_DIR/preprocess.py" \
        -r1 "$r1" -r2 "$r2" \
        -o "$output_path" -s "$sample_name" \
        -b1 "$whitelist_path" -b2 "$whitelist_path" \
        -l "$locus" -c "$cores" -q "$base_quality" \
        -bs "$preprocess_batch_size" \
        -cl "$compression_level" -cut "$cutadapt" \
        -l1 "$linker1" -l2 "$linker2" -m "$mm_rate" \
        -go "$gzip_output" \
        -cb "false" &> "$output_path/${sample_name}_preprocess.log" || {
            echo "Error: preprocessing failed for $sample_name; see $output_path/${sample_name}_preprocess.log" >&2
            exit 1
        }

    tmp_path=$(tail -n 1 $output_path/${sample_name}_preprocess.log)

    if [ "$preprocess_bc_ext" = "fq" ] && $gzip_after_enabled; then
        compress_fastq_file "$tmp_path/${sample_name}_bc_match_R1.fq" "$cores" "$compression_level" || exit 1
        compress_fastq_file "$tmp_path/${sample_name}_bc_match_R2.fq" "$cores" "$compression_level" || exit 1
    fi

    results=$output_path/results/$nonlocus_sample_name/$locus
    mkdir -p "$results"

    run_pixi python "$PYTHON_DIR/amplicon.py" \
        -bu "$tmp_path/${sample_name}_bc_match_R1.$bc_ext" \
        -dr "$tmp_path/${sample_name}_bc_match_R2.$bc_ext" \
        -o "$results" -d "$cutadapt" \
        --whitelist "$whitelist_path" \
        --sb-len "$sb_len" \
        --ub-len "$ub_len" \
        --umi_hd_threshold "$umi_hd_threshold" \
        --min-lb-len "$min_lb_len" \
        --initial-reads-cutoff "$initial_reads_cutoff" \
        --lb-error-rate "$lb_error_rate" \
        --lb-min-hd "$lb_min_hd" \
        --major-fraction-threshold-molecule "$major_fraction_threshold_molecule" \
        --reads-fraction-mode "$reads_fraction_mode" \
        --final-reads-cutoff "$reads_cutoff" \
        --slope-cutoff "$slope_cutoff" &> "$results/dbit.log" || {
            echo "Error: amplicon analysis failed for $sample_name; see $results/dbit.log" >&2
            exit 1
        }

    run_pixi python "$PYTHON_DIR/plot/heatmap.py" \
        -f "$results/final.csv" \
        -w "$whitelist_path" \
        -o "$results" || exit 1
done

if [ -n "$scratch" ]; then
    cp -r "$scratch_output"/* "$orig_output_path"/
    rm -rf "$scratch/amplicon"
fi
