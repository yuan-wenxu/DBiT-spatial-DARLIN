#!/bin/bash
set -o pipefail

SCRIPT_DIR=${STEP_SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)} || exit 1
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
for variable in x_spots_number y_spots_number length_spot interval \
    barcode_a_whitelist_path barcode_b_whitelist_path; do
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
    local status=$?
    trap - EXIT INT TERM HUP
    if [[ -n ${scratch_sample:-} && -d $scratch_sample ]]; then
        if (( status != 0 )) && [[ -n ${scratch_output:-} && -d $scratch_output && -n ${orig_output_path:-} ]]; then
            echo "Recovering mRNA scratch outputs after exit status $status: $orig_output_path" >&2
            mkdir -p "$orig_output_path" && cp -a "$scratch_output/." "$orig_output_path/" || \
                echo "Warning: failed to recover scratch outputs: $scratch_output" >&2
        fi
        rm -rf -- "$scratch_sample" || echo "Warning: failed to clean scratch directory: $scratch_sample" >&2
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

shopt -s nullglob
r1_files=("$fastq_path"/*_R1.fq.gz)
shopt -u nullglob
if (( ${#r1_files[@]} == 0 )); then
    echo "Error: mRNA FASTQ directory must contain at least one *_R1.fq.gz file: $fastq_path." >&2
    exit 1
fi
for fastq_file in "${r1_files[@]}"; do
    if [[ "$fastq_file" == *,* ]]; then
        echo "Error: mRNA FASTQ paths cannot contain commas because STAR uses commas to separate input files: $fastq_file" >&2
        exit 1
    fi
done
if [[ "$output_path" == *,* ]]; then
    echo "Error: mRNA output paths cannot contain commas because STAR uses commas to separate input files: $output_path" >&2
    exit 1
fi

mkdir -p "$output_path"
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
if [[ "$preprocess_bc_ext" == "fq" ]] && $gzip_after_enabled; then
    bc_ext="fq.gz"
fi

use_scratch=false
if [[ -n ${scratch:-} ]]; then
    use_scratch=true
    scratch_sample="$scratch/dbit/$run_id/mrna"
    scratch_input="$scratch_sample/input"
    scratch_output="$scratch_sample/output"
    enable_cleanup
fi

file_signature() {
    stat --printf='%n\t%s\t%Y\n' -- "$1"
}

preprocess_manifest_content() {
    local input_r1=$1
    local input_r2=$2
    printf 'manifest_version\t1\n'
    file_signature "$input_r1"
    file_signature "$input_r2"
    file_signature "$barcode_a_whitelist_path"
    file_signature "$barcode_b_whitelist_path"
    printf 'compression_level\t%s\n' "$compression_level"
    printf 'mm_rate\t%s\n' "$mm_rate"
    printf 'linker1\t%s\n' "$linker1"
    printf 'linker2\t%s\n' "$linker2"
    printf 'gzip_output\t%s\n' "$gzip_output"
    printf 'gzip_after_preprocess\t%s\n' "$gzip_after_preprocess"
    printf 'soloCBlen\t%s\n' "$soloCBlen"
    printf 'soloUMIlen\t%s\n' "$soloUMIlen"
}

manifest_matches() {
    local expected=$1
    local manifest_file=$2
    [[ -f "$manifest_file" ]] && cmp -s <(printf '%s\n' "$expected") "$manifest_file"
}

preprocessed_r1_files=()
preprocessed_r2_files=()

for r1 in "${r1_files[@]}"; do
    sample_name=$(basename "$r1" _R1.fq.gz)
    r2_orig="$fastq_path/${sample_name}_R2.fq.gz"
    if [[ ! -f "$r2_orig" ]]; then
        echo "Error: matching mRNA R2 file not found: $r2_orig" >&2
        exit 1
    fi

    log_file="$orig_output_path/${sample_name}_preprocess.log"
    pre_r1="$orig_output_path/fastq_umi_barcode/${sample_name}_bc_match_R1.$bc_ext"
    pre_r2="$orig_output_path/fastq_umi_barcode/${sample_name}_bc_match_R2.$bc_ext"
    preprocess_manifest="$orig_output_path/fastq_umi_barcode/${sample_name}_preprocess_manifest.tsv"
    current_preprocess_manifest=$(preprocess_manifest_content "$r1" "$r2_orig") || exit 1
    pre_done=false
    if [[ -s "$pre_r1" && -s "$pre_r2" && -f "$log_file" ]] && \
        manifest_matches "$current_preprocess_manifest" "$preprocess_manifest"; then
        pre_done=true
    fi

    # Step 1: preprocess each FASTQ pair independently.
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
            --reads1 "$step1_r1" --reads2 "$step1_r2" \
            --output "$step1_out" --sample "$sample_name" \
            --barcodeA_whitelist "$barcode_a_whitelist_path" \
            --barcodeB_whitelist "$barcode_b_whitelist_path" \
            --compression_level "$compression_level" \
            --mm_rate "$mm_rate" \
            --linker1 "$linker1" --linker2 "$linker2" \
            --correct_barcode "true" \
            --core "$preprocess_cores" \
            --batch_size "$preprocess_batch_size" \
            --gzip_output "$gzip_output" \
            --cb_len "$soloCBlen" \
            --umi_len "$soloUMIlen" 2>&1 | tee "$step1_log" || {
                echo "Error: preprocessing failed for $sample_name; see $step1_log" >&2
                exit 1
            }

        if [[ "$preprocess_bc_ext" == "fq" ]] && $gzip_after_enabled; then
            pre_dir="$step1_out/fastq_umi_barcode"
            compress_fastq_file "$pre_dir/${sample_name}_bc_match_R1.fq" "$preprocess_cores" "$compression_level" || exit 1
            compress_fastq_file "$pre_dir/${sample_name}_bc_match_R2.fq" "$preprocess_cores" "$compression_level" || exit 1
        fi

        # Keep step1 outputs in original output path for future skip checks.
        if $use_scratch; then
            cp "$step1_log" "$log_file"
            mkdir -p "$orig_output_path/fastq_umi_barcode"
            cp -a "$scratch_output/fastq_umi_barcode/." "$orig_output_path/fastq_umi_barcode/" || exit 1
        fi
        printf '%s\n' "$current_preprocess_manifest" > "$preprocess_manifest" || exit 1
    fi

    if [[ ! -s "$pre_r1" || ! -s "$pre_r2" ]]; then
        echo "Error: missing preprocess outputs for $sample_name." >&2
        exit 1
    fi
    preprocessed_r1_files+=("$pre_r1")
    preprocessed_r2_files+=("$pre_r2")
done

star_manifest_content() {
    printf 'manifest_version\t1\n'
    printf 'genome_dir\t%s\n' "$genome_dir"
    printf 'soloCBstart\t%s\n' "$soloCBstart"
    printf 'soloCBlen\t%s\n' "$soloCBlen"
    printf 'soloUMIstart\t%s\n' "$soloUMIstart"
    printf 'soloUMIlen\t%s\n' "$soloUMIlen"
    local index
    for index in "${!preprocessed_r1_files[@]}"; do
        printf 'pair\t%s\n' "$index"
        preprocess_manifest_content \
            "${r1_files[index]}" \
            "${r1_files[index]%_R1.fq.gz}_R2.fq.gz"
        file_signature "${preprocessed_r1_files[index]}"
        file_signature "${preprocessed_r2_files[index]}"
    done
}

current_star_manifest=$(star_manifest_content) || exit 1
star_manifest="$final_results/star_input_manifest.tsv"
star_done=false
if star_outputs_complete "$final_results" && \
    manifest_matches "$current_star_manifest" "$star_manifest"; then
    star_done=true
elif [[ -d "$final_results" ]]; then
    echo "Warning: incomplete STAR outputs or changed mRNA inputs found; rerunning Step2." >&2
fi

# Step 2: combine all preprocessed pairs in one STAR run.
if $star_done; then
    echo "Step2 STAR already done for ${#r1_files[@]} FASTQ pair(s), skipping..."
else
    if $use_scratch; then
        star_input="$scratch_sample/star_input"
        mkdir -p "$star_input" "$scratch_output"
        star_r1_files=()
        star_r2_files=()
        for index in "${!preprocessed_r1_files[@]}"; do
            cp "${preprocessed_r1_files[index]}" "${preprocessed_r2_files[index]}" "$star_input/" || exit 1
            star_r1_files+=("$star_input/$(basename "${preprocessed_r1_files[index]}")")
            star_r2_files+=("$star_input/$(basename "${preprocessed_r2_files[index]}")")
        done
        star_results="$scratch_output/results"
        # Only clean our own run_id-scoped scratch directory (safe for concurrency).
        if [[ -d "$star_results" ]]; then
            echo "Removing stale scratch STAR outputs: $star_results"
            rm -rf -- "$star_results" || {
                echo "Error: failed to remove stale STAR outputs: $star_results" >&2
                exit 1
            }
        fi
    else
        star_r1_files=("${preprocessed_r1_files[@]}")
        star_r2_files=("${preprocessed_r2_files[@]}")
        star_results="$final_results"
        # Non-scratch mode: clean shared final_results (not concurrency-safe; use scratch for parallel runs).
        if [[ -d "$final_results" ]]; then
            echo "Removing incomplete or outdated STAR outputs while preserving deconv: $final_results"
            find "$final_results" -mindepth 1 -maxdepth 1 ! -name deconv \
                -exec rm -rf -- {} + || {
                echo "Error: failed to remove incomplete STAR outputs: $final_results" >&2
                exit 1
            }
        fi
    fi
    mkdir -p "$star_results"
    star_read_args=()
    if [[ "$bc_ext" == "fq.gz" ]]; then
        star_read_args=(--readFilesCommand zcat)
    fi

    star_r1_csv=$(IFS=,; printf '%s' "${star_r1_files[*]}")
    star_r2_csv=$(IFS=,; printf '%s' "${star_r2_files[*]}")
    echo "Running STAR with ${#star_r1_files[@]} FASTQ pair(s)."
    run_pixi STAR \
        --runMode alignReads \
        --runThreadN "$star_threads" \
        --genomeDir "$genome_dir" \
        "${star_read_args[@]}" \
        --readFilesIn "$star_r2_csv" "$star_r1_csv" \
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
            echo "Error: STAR failed; see $star_results/STAR.log" >&2
            exit 1
        }

    if ! star_outputs_complete "$star_results"; then
        echo "Error: STAR finished without all required outputs: $star_results" >&2
        exit 1
    fi
    printf '%s\n' "$current_star_manifest" > "$star_results/star_input_manifest.tsv" || exit 1

    # Step3 is local only, so copy step2 results back first when using scratch.
    if $use_scratch; then
        mkdir -p "$orig_output_path/results"
        cp -a "$star_results/." "$orig_output_path/results/" || {
            echo "Error: failed to copy STAR results from scratch" >&2
            exit 1
        }
    fi
fi

# Step 3: always run locally, no skip and no scratch.
if ! star_outputs_complete "$final_results"; then
    echo "Error: incomplete STAR outputs for mRNA QC: $final_results" >&2
    exit 1
fi
run_pixi python "$PYTHON_DIR/mrna.py" \
    --file_path "$final_results/Solo.out" \
    --barcodeA_whitelist "$barcode_a_whitelist_path" \
    --barcodeB_whitelist "$barcode_b_whitelist_path" \
    --umi_min "$umi_min" --gene_min "$gene_min" --min_cells "$min_cells" \
    --cb-len "$soloCBlen" \
    --x_spots_number "$x_spots_number" --y_spots_number "$y_spots_number" \
    --length_spot "$length_spot" --interval "$interval" \
    --pixel_length "$pixel_length" 2>&1 | tee "$final_results/Solo.out/qc.log" || {
        echo "Error: mRNA QC failed; see $final_results/Solo.out/qc.log" >&2
        exit 1
    }

if $use_scratch; then
    rm -rf -- "$scratch_sample"
    scratch_sample=""
fi
