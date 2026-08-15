#!/bin/bash
set -o pipefail

SCRIPT_DIR=${QC_SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)} || exit 1
REPO_DIR=${REPO_DIR:-$(cd "$SCRIPT_DIR/../.." && pwd)} || exit 1
PYTHON_DIR="$SCRIPT_DIR/python"

if [[ $# -ne 1 ]]; then echo "Error: expected one config file argument" >&2; exit 1; fi
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

run_id=${SLURM_JOB_ID:-amplicon_$$}
scratch_run_dir=""

cleanup_scratch() {
    local status=$?
    trap - EXIT INT TERM HUP
    if [[ -n ${scratch_run_dir:-} && -d $scratch_run_dir ]]; then
        if (( status != 0 )) && [[ -n ${scratch_output:-} && -d $scratch_output && -n ${orig_output_path:-} ]]; then
            echo "Recovering amplicon scratch outputs after exit status $status: $orig_output_path" >&2
            mkdir -p "$orig_output_path" && cp -a "$scratch_output/." "$orig_output_path/" || \
                echo "Warning: failed to recover scratch outputs: $scratch_output" >&2
        fi
        rm -rf -- "$scratch_run_dir" || echo "Warning: failed to clean scratch directory: $scratch_run_dir" >&2
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
        pixi run -e "$pixi_env" "$@"
    )
}

compress_fastq_file() {
    local fq="$1"
    local threads="$2"
    local level="$3"
    run_pixi pigz -f -p "$threads" "-$level" "$fq"
}

run_amplicon_cutadapt() {
    local locus="$1"
    local reads1="$2"
    local reads2="$3"
    local sample="$4"
    local trim_dir="$output_path/fastq_cut"
    local prime3 prime5

    case "$locus" in
        CA)
            prime3="AGAATTCTAACTAGA"
            prime5="GTACAAGTAAGCGGC"
            ;;
        RA)
            prime3="GTCTGCTGTGTGCCT"
            prime5="ACAAGTAAAGCGGCC"
            ;;
        TA)
            prime3="GTCTTGTCGGTGCCT"
            prime5="TCGGTACCTCGCGAA"
            ;;
        *)
            echo "Error: unsupported DARLIN locus for cutadapt: $locus" >&2
            return 1
            ;;
    esac

    mkdir -p "$trim_dir"
    local temporary_r1="$trim_dir/${sample}_R1.tmp.fq.gz"
    local temporary_r2="$trim_dir/${sample}_R2.tmp.fq.gz"
    preprocess_r1="$trim_dir/${sample}_R1.trimmed.fq.gz"
    preprocess_r2="$trim_dir/${sample}_R2.trimmed.fq.gz"

    run_pixi cutadapt -j "$cores" -A "$prime3" \
        --overlap "${#prime3}" --error-rate 0 \
        -q "$base_quality" -m 10 --discard-untrimmed \
        -o "$temporary_r1" -p "$temporary_r2" "$reads1" "$reads2" || return 1
    run_pixi cutadapt -j "$cores" -G "$prime5" \
        --overlap "${#prime5}" --error-rate 0 \
        -q "$base_quality" -m 10 --discard-untrimmed \
        -o "$preprocess_r1" -p "$preprocess_r2" \
        "$temporary_r1" "$temporary_r2" || return 1
    rm -f -- "$temporary_r1" "$temporary_r2"
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

case "${cutadapt,,}" in
    true|yes|1) cutadapt_enabled=true ;;
    false|no|0) cutadapt_enabled=false ;;
    *)
        echo "Error: cutadapt must be true or false" >&2
        exit 1
        ;;
esac

if [ -n "$scratch" ]; then
    scratch_input="$scratch/dbit/$run_id/amplicon/input"
    scratch_output="$scratch/dbit/$run_id/amplicon/output"
    scratch_run_dir="$scratch/dbit/$run_id/amplicon"
    enable_cleanup
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

    # Cutadapt and extract UMI and barcode
    preprocess_r1="$r1"
    preprocess_r2="$r2"
    if $cutadapt_enabled; then
        run_amplicon_cutadapt \
            "$locus" "$r1" "$r2" "$sample_name" || {
                echo "Error: cutadapt failed for $sample_name" >&2
                exit 1
            }
    fi

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
        -r1 "$preprocess_r1" -r2 "$preprocess_r2" \
        -o "$output_path" -s "$sample_name" \
        -b1 "$whitelist_path" -b2 "$whitelist_path" \
        -c "$cores" \
        -bs "$preprocess_batch_size" \
        -cl "$compression_level" \
        -l1 "$linker1" -l2 "$linker2" -m "$mm_rate" \
        -go "$gzip_output" \
        -cb "false" \
        -bl "$sb_len" \
        -ul "$ub_len" 2>&1 | tee "$output_path/${sample_name}_preprocess.log" || {
            echo "Error: preprocessing failed for $sample_name; see $output_path/${sample_name}_preprocess.log" >&2
            exit 1
        }

    tmp_path=$(tail -n 1 $output_path/${sample_name}_preprocess.log)

    if [ "$preprocess_bc_ext" = "fq" ] && $gzip_after_enabled; then
        compress_fastq_file "$tmp_path/${sample_name}_bc_match_R1.fq" "$cores" "$compression_level" || exit 1
        compress_fastq_file "$tmp_path/${sample_name}_bc_match_R2.fq" "$cores" "$compression_level" || exit 1
    fi

    results="$output_path/results/$locus"
    mkdir -p "$results"

    run_pixi python "$PYTHON_DIR/amplicon.py" \
        -bu "$tmp_path/${sample_name}_bc_match_R1.$bc_ext" \
        -dr "$tmp_path/${sample_name}_bc_match_R2.$bc_ext" \
        -o "$results" -d "$cutadapt" \
        --whitelist "$whitelist_path" \
        --sb-len "$sb_len" \
        --ub-len "$ub_len" \
        --x-spots-number "$x_spots_number" \
        --y-spots-number "$y_spots_number" \
        --umi_hd_threshold "$umi_hd_threshold" \
        --min-lb-len "$min_lb_len" \
        --initial-reads-cutoff "$initial_reads_cutoff" \
        --lb-error-rate "$lb_error_rate" \
        --lb-min-hd "$lb_min_hd" \
        --major-fraction-threshold-molecule "$major_fraction_threshold_molecule" \
        --reads-fraction-mode "$reads_fraction_mode" \
        --final-reads-cutoff "$reads_cutoff" \
        --slope-cutoff "$slope_cutoff" 2>&1 | tee "$results/dbit.log" || {
            echo "Error: amplicon analysis failed for $sample_name; see $results/dbit.log" >&2
            exit 1
        }

done

if [ -n "$scratch" ]; then
    mkdir -p "$orig_output_path"
    cp -a "$scratch_output/." "$orig_output_path/"
fi
