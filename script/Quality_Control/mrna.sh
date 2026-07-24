#!/bin/bash
set -o pipefail

show_help() {
    cat << EOF
Usage: $0 <config_file>

Process mRNA sequencing data with preprocessing, STAR alignment, and quality control.

Arguments:
  config_file   Per-dataset QC configuration file

Examples:
  $0 dbit.config.sh
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
fastq_path=$mrna_fastq_path
output_path=${mrna_output_path}
preprocess_cores=${mrna_cores}
umi_min=${umi_min:-900}
gene_min=${gene_min:-300}
min_cells=${min_cells:-3}

if [[ -z ${mrna_fastq_path:-} || -z ${genome_dir:-} ]]; then
    echo "Error: mrna_fastq_path and genome_dir must be set in the QC config." >&2
    exit 1
fi
for variable in x_spots_number y_spots_number length_spot interval whitelist_path; do
    if [[ -z ${!variable:-} ]]; then
        echo "Run this script through dbit.sh so --chip is resolved." >&2
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

run_id=${SLURM_JOB_ID:-mrna_$$}
scratch_sample=""

cleanup_scratch() {
    local run_dir="${scratch:-}/dbit/${run_id:-}"
    if [[ -n "${run_id:-}" && -n "${scratch:-}" && -d "$run_dir" ]]; then
        rm -rf -- "$run_dir"
    fi
}

trap cleanup_scratch EXIT
trap 'exit 130' INT
trap 'exit 143' TERM HUP

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

star_outputs_complete() {
    local result_dir=$1
    local required_file
    local required_files=(
        "Aligned.sortedByCoord.out.bam"
        "Log.final.out"
        "Solo.out/GeneFull/raw/matrix.mtx"
        "Solo.out/GeneFull/raw/barcodes.tsv"
        "Solo.out/GeneFull/raw/features.tsv"
    )
    for required_file in "${required_files[@]}"; do
        if [[ ! -s "$result_dir/$required_file" ]]; then
            return 1
        fi
    done
    return 0
}

# Validate inputs
fastq_path=$(normalize_dir_path "$fastq_path")
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

orig_output_path="$output_path"

mkdir -p "$output_path"

for r1 in "$fastq_path"/*_R1.fq.gz; do
    [ -e "$r1" ] || { echo "Error: no *_R1.fq.gz files found in $fastq_path" >&2; exit 1; }
    sample_name=$(basename "$r1" | sed 's/_R1.fq.gz//')
    r2_orig="$fastq_path/${sample_name}_R2.fq.gz"

    log_file="$orig_output_path/${sample_name}_preprocess.log"
    final_results="$orig_output_path/results"

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

    pre_file="$orig_output_path/fastq_umi_barcode/${sample_name}_bc_match_R1.$bc_ext"
    pre_done=false
    star_done=false
    if [ -f "$pre_file" ] && [ -f "$log_file" ]; then pre_done=true; fi
    if star_outputs_complete "$final_results"; then
        star_done=true
    elif [[ -d "$final_results" ]]; then
        echo "Warning: incomplete STAR outputs found for $sample_name; rerunning Step2." >&2
    fi

    use_scratch=false
    if [ -n "$scratch" ]; then
        use_scratch=true
        scratch_sample="$scratch/dbit/$run_id/mrna"
        scratch_input="$scratch_sample/input"
        scratch_output="$scratch_sample/output"
    fi

    # Step 1: preprocess (skip if already done)
    if $pre_done; then
        echo "Step1 preprocess already done for $sample_name, skipping..."
    else
        if $use_scratch; then
            mkdir -p "$scratch_input" "$scratch_output"
            cp "$r1" "$r2_orig" "$scratch_input/"
            step1_r1="$scratch_input/${sample_name}_R1.fq.gz"
            step1_r2="$scratch_input/${sample_name}_R2.fq.gz"
            step1_out="$scratch_output"
            step1_log="$scratch_output/${sample_name}_preprocess.log"
        else
            step1_r1="$r1"
            step1_r2="$r2_orig"
            step1_out="$orig_output_path"
            step1_log="$log_file"
        fi

        run_pixi python "$PYTHON_DIR/preprocess.py" \
            -r1 "$step1_r1" -r2 "$step1_r2" \
            -o "$step1_out" -s "$sample_name" \
            -b1 "$whitelist_path" -b2 "$whitelist_path" \
            -cl "$compression_level" -m "$mm_rate" \
            -l1 "$linker1" -l2 "$linker2" -cb "true" \
            -c "$preprocess_cores" \
            -bs "$preprocess_batch_size" \
            -go "$gzip_output" \
            -bmd "$bc_max_dist" 2>&1 | tee "$step1_log" || {
                echo "Error: preprocessing failed for $sample_name; see $step1_log" >&2
                exit 1
            }

        if [ "$preprocess_bc_ext" = "fq" ] && $gzip_after_enabled; then
            pre_dir="$step1_out/fastq_umi_barcode"
            compress_fastq_file "$pre_dir/${sample_name}_bc_match_R1.fq" "$preprocess_cores" "$compression_level" || exit 1
            compress_fastq_file "$pre_dir/${sample_name}_bc_match_R2.fq" "$preprocess_cores" "$compression_level" || exit 1
        fi

        # Keep step1 outputs in original output path for future skip checks.
        if $use_scratch; then
            cp "$step1_log" "$log_file"
            cp -r "$scratch_output"/fastq_umi_barcode "$orig_output_path"/ 2>/dev/null || true
        fi
    fi

    # Resolve step2 input from existing/preprocessed files in original output path.
    pre_file="$orig_output_path/fastq_umi_barcode/${sample_name}_bc_match_R1.$bc_ext"
    if [ -f "$pre_file" ]; then
        pre_r1="$pre_file"
        pre_r2="${pre_r1%_R1.$bc_ext}_R2.$bc_ext"
        tmp_path="$(dirname "$pre_r1")"
    else
        echo "Error: missing preprocess outputs for $sample_name." >&2
        exit 1
    fi

    # Step 2: STAR (skip if already done)
    if $star_done; then
        echo "Step2 STAR already done for $sample_name, skipping..."
    else
        if $use_scratch; then
            mkdir -p "$scratch_input" "$scratch_output"
            cp "$pre_r1" "$pre_r2" "$scratch_input/"
            star_input="$scratch_input"
            star_results="$scratch_output/results"
            # Only clean our own run_id-scoped scratch directory (safe for concurrency)
            if [[ -d "$star_results" ]]; then
                echo "Removing stale scratch STAR outputs for $sample_name: $star_results"
                rm -rf -- "$star_results" || {
                    echo "Error: failed to remove stale STAR outputs: $star_results" >&2
                    exit 1
                }
            fi
        else
            star_input="$tmp_path"
            star_results="$orig_output_path/results"
            # Non-scratch mode: clean shared final_results (not concurrency-safe; use scratch for parallel runs)
            if [[ -d "$final_results" ]]; then
                echo "Removing incomplete STAR outputs for $sample_name while preserving deconv: $final_results"
                find "$final_results" -mindepth 1 -maxdepth 1 ! -name deconv \
                    -exec rm -rf -- {} + || {
                    echo "Error: failed to remove incomplete STAR outputs: $final_results" >&2
                    exit 1
                }
            fi
        fi
        mkdir -p "$star_results"
        star_read_args=()
        if [ "$bc_ext" = "fq.gz" ]; then
            star_read_args=(--readFilesCommand zcat)
        fi

        run_pixi STAR \
            --runMode alignReads \
            --runThreadN "$star_threads" \
            --genomeDir "$genome_dir" \
            "${star_read_args[@]}" \
            --readFilesIn "$star_input/${sample_name}_bc_match_R2.$bc_ext" "$star_input/${sample_name}_bc_match_R1.$bc_ext" \
            --outFileNamePrefix "$star_results/" \
            --outTmpDir "$star_results/solotmp" \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS nM CB UB GX GN \
            --soloType CB_UMI_Simple \
            --soloCBstart "$soloCBstart" \
            --soloCBlen "$soloCBlen" \
            --soloUMIstart "$soloUMIstart" \
            --soloUMIlen "$soloUMIlen" \
            --soloBarcodeReadLength 0 \
            --soloCBwhitelist None \
            --soloCellFilter None \
            --soloFeatures GeneFull \
            --bamRemoveDuplicatesType UniqueIdentical \
            --quantMode GeneCounts 2>&1 | tee "$star_results/STAR.log" || {
                echo "Error: STAR failed for $sample_name; see $star_results/STAR.log" >&2
                exit 1
            }

        if ! star_outputs_complete "$star_results"; then
            echo "Error: STAR finished without all required outputs for $sample_name: $star_results" >&2
            exit 1
        fi

        # Step3 is local only, so copy step2 results back first when using scratch.
        if $use_scratch; then
            mkdir -p "$orig_output_path/results"
            cp -r "$star_results"/* "$orig_output_path/results"/ || {
                echo "Error: failed to copy STAR results from scratch for $sample_name" >&2
                exit 1
            }
        fi
    fi

    # Step 3: always run locally, no skip and no scratch.
    if ! star_outputs_complete "$final_results"; then
        echo "Error: incomplete STAR outputs for mRNA QC: $final_results" >&2
        exit 1
    fi
    run_pixi python "$PYTHON_DIR/mrna.py" -f "$final_results/Solo.out" -w "$whitelist_path" \
        -umi_min "$umi_min" -gene_min "$gene_min" -min_cells "$min_cells" \
        --x_spots_number "$x_spots_number" --y_spots_number "$y_spots_number" \
        --length_spot "$length_spot" --interval "$interval" \
        --pixel_length "$pixel_length" 2>&1 | tee "$final_results/Solo.out/qc.log" || {
            echo "Error: mRNA QC failed for $sample_name; see $final_results/Solo.out/qc.log" >&2
            exit 1
        }

    if $use_scratch; then
        rm -rf -- "$scratch_sample"
        scratch_sample=""
    fi
done
