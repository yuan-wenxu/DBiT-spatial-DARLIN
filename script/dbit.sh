#!/bin/bash
set -o pipefail

SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || exit 1
SCRIPT_DIR=$(cd "$(dirname "$SCRIPT_PATH")" && pwd) || exit 1
REPO_DIR=$(cd "$SCRIPT_DIR/.." && pwd) || exit 1
START_DIR=$(pwd -P) || exit 1
STEP_SCRIPT_DIR="$SCRIPT_DIR/step"
PROGRAM_NAME=$(basename "$0")
CHIP_FILE="$REPO_DIR/config/chip.sh"
export STEP_SCRIPT_DIR REPO_DIR

if [[ ! -f "$CHIP_FILE" ]]; then
    echo "Error: chip preset file not found: $CHIP_FILE" >&2
    exit 1
fi
# shellcheck source=/dev/null
source "$CHIP_FILE"

show_help() {
    cat <<EOF
Usage: $PROGRAM_NAME <step> [options]

Steps:
  init          Initialize or display configuration in the current directory
  mrna          Process transcriptome FASTQs and run spatial mRNA QC
  saturation    Downsample mRNA FASTQs and run mRNA QC at each fraction
  darlin        Process DARLIN FASTQs
  image         Create a tissue mask and apply tissue filtering

Run '$PROGRAM_NAME <step> -h' to show parameters for one step.
EOF
}

show_init_help() {
    cat <<EOF
Usage: $PROGRAM_NAME init

Initialize dbit.config.sh in the current directory.
If dbit.config.sh already exists, display its contents instead.
EOF
}

show_mrna_help() {
    cat <<EOF
Usage: $PROGRAM_NAME mrna [--config <file>] [options]

Optional:
  --config <file>        Configuration file (default: ./dbit.config.sh)
  --input <path>         Transcriptome FASTQ directory; required only before stored
  --chip <name>          $(chip_preset_names_csv); required only before stored
  --umi-min <int>        Non-negative minimum UMI count per spot (default: 900)
  --gene-min <int>       Non-negative minimum gene count per spot (default: 300)
  --min-cell <int>       Positive minimum cells per gene (default: 3)
EOF
}

show_darlin_help() {
    cat <<EOF
Usage: $PROGRAM_NAME darlin [--config <file>] [options]

Optional:
  --config <file>                                 Configuration file (default: ./dbit.config.sh)
  --input <path>                                  DARLIN FASTQ directory; required only before stored
  --chip <name>                                   $(chip_preset_names_csv); required only before stored
  --initial-reads-cutoff <int>                    Non-negative reads cutoff for initial filtering (default: 100)
  --major-fraction-threshold-molecule <float>     Major-molecule reads fraction from 0 to 1 (default: 0.8)
  --reads-fraction-mode <sum|max>                 Mode for calculating reads fraction (default: sum)
  --reads-cutoff <int>                            Non-negative reads cutoff for final filtering (default: 10)
  --slope-cutoff <float>                          Non-negative slope cutoff for final filtering (default: 10)
EOF
}

show_saturation_help() {
    cat <<EOF
Usage: $PROGRAM_NAME saturation [--config <file>] [options]

The input is read from mrna_fastq_path stored by the mrna step. Downsampled
FASTQs and mRNA results are written below <mRNA FASTQ parent>/saturation/.

Optional:
  --config <file>       Configuration file (default: ./dbit.config.sh)
  --fractions <list>    Comma-separated fractions (default: 0.01,0.02,0.05,0.1,0.2,0.5)
EOF
}

show_image_help() {
    cat <<EOF
Usage: $PROGRAM_NAME image [--config <file>] [options]

Optional:
  --config <file>         Configuration file (default: ./dbit.config.sh)
  --input <path>          Full-resolution image; adjacent mask.png defines the frame
  --chip <name>           $(chip_preset_names_csv); required only before stored
EOF
}

show_step_help() {
    case "$1" in
        init) show_init_help ;;
        mrna) show_mrna_help ;;
        saturation) show_saturation_help ;;
        darlin) show_darlin_help ;;
        image) show_image_help ;;
    esac
}

if [[ $# -eq 0 || ${1:-} == -h || ${1:-} == --help ]]; then show_help; exit 0; fi

step=$1
if [[ ${2:-} == -h || ${2:-} == --help ]]; then
    case "$step" in
        init|mrna|saturation|darlin|image) show_step_help "$step"; exit 0 ;;
    esac
fi
shift

case "$step" in
    init|mrna|saturation|darlin|image) ;;
    *)
        echo "Error: unsupported step '$step'." >&2
        echo "Valid steps: init, mrna, saturation, darlin, image." >&2
        exit 1
        ;;
esac

require_option_value() {
    if [[ $# -lt 2 || $2 == --* ]]; then
        echo "Error: option '$1' requires a value." >&2
        exit 1
    fi
}

# Handle init step early (no config required)
if [[ "$step" == "init" ]]; then
    [[ ${2:-} == -h || ${2:-} == --help ]] && { show_init_help; exit 0; }

    init_config_file="$START_DIR/dbit.config.sh"

    if [[ -f "$init_config_file" ]]; then
        echo "Config file already exists: $(realpath "$init_config_file")"
        echo "---"
        cat "$init_config_file"
    else
        if [[ -f "$REPO_DIR/config/dbit.config.sh" ]]; then
            cp "$REPO_DIR/config/dbit.config.sh" "$init_config_file"
        else
            cp "$REPO_DIR/config/dbit.config.example.sh" "$init_config_file"
        fi
        echo "Created config file: $(realpath "$init_config_file")"
    fi
    exit 0
fi

input_path=""
input_from_cli=false
config_file=""
selected_chip=""
chip_from_cli=false
cli_umi_min=""
cli_gene_min=""
cli_min_cell=""
cli_initial_reads_cutoff=""
cli_major_fraction_threshold_molecule=""
cli_reads_fraction_mode=""
cli_reads_cutoff=""
cli_slope_cutoff=""
cli_saturation_fractions=""

require_step_option() {
    local option=$1
    shift
    local allowed_step
    for allowed_step in "$@"; do
        [[ "$step" == "$allowed_step" ]] && return 0
    done
    echo "Error: option '$option' is not valid for the '$step' step." >&2
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config) require_step_option "$1" mrna saturation darlin image; require_option_value "$@"; config_file=$2; shift 2 ;;
        --input) require_step_option "$1" mrna darlin image; require_option_value "$@"; input_path=$2; input_from_cli=true; shift 2 ;;
        --chip) require_step_option "$1" mrna darlin image; require_option_value "$@"; selected_chip=$2; chip_from_cli=true; shift 2 ;;
        --umi-min) require_step_option "$1" mrna; require_option_value "$@"; cli_umi_min=$2; shift 2 ;;
        --gene-min) require_step_option "$1" mrna; require_option_value "$@"; cli_gene_min=$2; shift 2 ;;
        --min-cell) require_step_option "$1" mrna; require_option_value "$@"; cli_min_cell=$2; shift 2 ;;
        --fractions) require_step_option "$1" saturation; require_option_value "$@"; cli_saturation_fractions=$2; shift 2 ;;
        --initial-reads-cutoff) require_step_option "$1" darlin; require_option_value "$@"; cli_initial_reads_cutoff=$2; shift 2 ;;
        --major-fraction-threshold-molecule) require_step_option "$1" darlin; require_option_value "$@"; cli_major_fraction_threshold_molecule=$2; shift 2 ;;
        --reads-fraction-mode) require_step_option "$1" darlin; require_option_value "$@"; cli_reads_fraction_mode=$2; shift 2 ;;
        --reads-cutoff) require_step_option "$1" darlin; require_option_value "$@"; cli_reads_cutoff=$2; shift 2 ;;
        --slope-cutoff) require_step_option "$1" darlin; require_option_value "$@"; cli_slope_cutoff=$2; shift 2 ;;
        -h|--help) show_step_help "$step"; exit 0 ;;
        *) echo "Error: unknown option or argument '$1'." >&2; exit 1 ;;
    esac
done

validate_nonnegative_integer() {
    local option=$1
    local value=$2
    if [[ -n "$value" && ! "$value" =~ ^[0-9]+$ ]]; then
        echo "Error: $option must be a non-negative integer; got '$value'." >&2
        exit 1
    fi
}

validate_positive_integer() {
    local option=$1
    local value=$2
    if [[ -n "$value" && ( ! "$value" =~ ^[0-9]+$ || "$value" =~ ^0+$ ) ]]; then
        echo "Error: $option must be a positive integer; got '$value'." >&2
        exit 1
    fi
}

validate_nonnegative_number() {
    local option=$1
    local value=$2
    if [[ -n "$value" && ! "$value" =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$ ]]; then
        echo "Error: $option must be a non-negative number; got '$value'." >&2
        exit 1
    fi
}

validate_positive_number() {
    local option=$1
    local value=$2
    if [[ -n "$value" && ! "$value" =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$ ]]; then
        echo "Error: $option must be a positive number; got '$value'." >&2
        exit 1
    fi
    if [[ -n "$value" ]] && ! awk -v value="$value" 'BEGIN { exit !(value > 0) }'; then
        echo "Error: $option must be greater than zero; got '$value'." >&2
        exit 1
    fi
}

validate_fraction() {
    local option=$1
    local value=$2
    if [[ -n "$value" && ! "$value" =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$ ]]; then
        echo "Error: $option must be a number from 0 to 1; got '$value'." >&2
        exit 1
    fi
    if [[ -n "$value" ]] && ! awk -v value="$value" 'BEGIN { exit !(value >= 0 && value <= 1) }'; then
        echo "Error: $option must be a number from 0 to 1; got '$value'." >&2
        exit 1
    fi
}

validate_single_mrna_fastq_pair() {
    local input_dir=$1
    local -a r1_files=()
    local r2_file

    mapfile -d '' -t r1_files < <(
        find "$input_dir" -maxdepth 1 -type f -name '*_R1.fq.gz' -print0
    )
    if (( ${#r1_files[@]} != 1 )); then
        echo "Error: mRNA FASTQ directory must contain exactly one *_R1.fq.gz file; found ${#r1_files[@]} in $input_dir." >&2
        exit 1
    fi

    r2_file="${r1_files[0]%_R1.fq.gz}_R2.fq.gz"
    if [[ ! -f "$r2_file" ]]; then
        echo "Error: matching mRNA R2 file not found: $r2_file" >&2
        exit 1
    fi
}

validate_nonnegative_integer --umi-min "$cli_umi_min"
validate_nonnegative_integer --gene-min "$cli_gene_min"
validate_positive_integer --min-cell "$cli_min_cell"
validate_nonnegative_integer --initial-reads-cutoff "$cli_initial_reads_cutoff"
validate_fraction --major-fraction-threshold-molecule "$cli_major_fraction_threshold_molecule"
validate_nonnegative_integer --reads-cutoff "$cli_reads_cutoff"
validate_nonnegative_number --slope-cutoff "$cli_slope_cutoff"
if [[ -n "$cli_saturation_fractions" ]]; then
    IFS=',' read -ra saturation_fraction_values <<< "$cli_saturation_fractions"
    if [[ ${#saturation_fraction_values[@]} -eq 0 ]]; then
        echo "Error: --fractions requires at least one fraction." >&2
        exit 1
    fi
    for saturation_fraction in "${saturation_fraction_values[@]}"; do
        saturation_fraction=${saturation_fraction//[[:space:]]/}
        validate_positive_number --fractions "$saturation_fraction"
        if ! awk -v value="$saturation_fraction" 'BEGIN { exit !(value <= 1) }'; then
            echo "Error: each --fractions value must be at most 1; got '$saturation_fraction'." >&2
            exit 1
        fi
    done
fi
if [[ -n "$cli_reads_fraction_mode" ]]; then
    case "$cli_reads_fraction_mode" in
        sum|max) ;;
        *)
            echo "Error: --reads-fraction-mode must be sum or max; got '$cli_reads_fraction_mode'." >&2
            exit 1
            ;;
    esac
fi

if [[ -z "$config_file" ]]; then
    config_file="$START_DIR/dbit.config.sh"
fi
if [[ ! -f "$config_file" ]]; then
    echo "Error: config file not found: $config_file" >&2
    echo "Run 'dbit init' to create dbit.config.sh, or pass --config <file>." >&2
    exit 1
fi
config_abs=$(realpath "$config_file")
echo "Using config: $config_abs"

# Load stored values first. Command-line values below take precedence and are
# appended to the config so that later runs can reuse them.
source "$config_abs"

set_config_value() {
    local key=$1
    local value=$2
    printf '%s=%q\n' "$key" "$value" >> "$config_abs"
}

if ! $chip_from_cli; then selected_chip=${chip:-}; fi
if [[ -n "$selected_chip" ]] && ! chip_preset_is_supported "$selected_chip"; then
    echo "Error: unsupported chip '$selected_chip'." >&2
    echo "Valid chips: $(chip_preset_names_csv)." >&2
    exit 1
fi

stored_input=""
case "$step" in
    mrna) stored_input=${mrna_fastq_path:-} ;;
    saturation) stored_input=${mrna_fastq_path:-} ;;
    darlin) stored_input=${darlin_fastq_path:-} ;;
    image) stored_input=${image_path:-} ;;
esac
[[ -n "$input_path" ]] || input_path=$stored_input
if [[ "$step" =~ ^(mrna|saturation|darlin|image)$ && -z "$input_path" ]]; then
    if [[ "$step" == saturation ]]; then
        echo "Error: mrna_fastq_path is not stored in $config_abs; run 'dbit mrna --input <fastq_dir>' first." >&2
    else
        echo "Error: --input is required the first time; no input path for '$step' is stored in $config_abs." >&2
    fi
    exit 1
fi
if [[ "$step" =~ ^(mrna|darlin|image)$ && -z "$selected_chip" ]]; then
    echo "Error: --chip is required the first time; no chip is stored in $config_abs." >&2
    echo "Run '$PROGRAM_NAME $step' again with --chip <name>. Valid chips: $(chip_preset_names_csv)." >&2
    exit 1
fi
if [[ "$step" == mrna ]]; then
    if [[ -z ${genome_dir:-} ]]; then
        echo "Error: genome_dir must be set in $config_abs." >&2
        exit 1
    fi
fi
if [[ "$step" == mrna || "$step" == saturation || "$step" == darlin || "$step" == image ]]; then
    input_abs=$(realpath -m "$input_path")
fi

case "$step" in
    mrna)
        validate_single_mrna_fastq_pair "$input_abs"
        output_path=$(dirname "$input_abs")
        if $input_from_cli; then
            set_config_value mrna_fastq_path "$input_abs"
            set_config_value mrna_output_path "$output_path"
            set_config_value mrna_dir "$output_path/results/Solo.out/GeneFull"
        fi
        [[ -n "$cli_umi_min" ]] && set_config_value umi_min "$cli_umi_min"
        [[ -n "$cli_gene_min" ]] && set_config_value gene_min "$cli_gene_min"
        [[ -n "$cli_min_cell" ]] && set_config_value min_cells "$cli_min_cell"
        ;;
    saturation)
        validate_single_mrna_fastq_pair "$input_abs"
        ;;
    darlin)
        first_r1=$(find "$input_abs" -maxdepth 1 -type f -name '*_R1.fq.gz' -print -quit)
        if [[ -z "$first_r1" ]]; then
            echo "Error: no *_R1.fq.gz file found in $input_abs" >&2
            exit 1
        fi
        output_path=$(dirname "$input_abs")
        if $input_from_cli; then
            set_config_value darlin_fastq_path "$input_abs"
            set_config_value darlin_output_path "$output_path"
            set_config_value darlin_dir "$output_path/results"
        fi
        [[ -n "$cli_initial_reads_cutoff" ]] && set_config_value initial_reads_cutoff "$cli_initial_reads_cutoff"
        [[ -n "$cli_major_fraction_threshold_molecule" ]] && set_config_value major_fraction_threshold_molecule "$cli_major_fraction_threshold_molecule"
        [[ -n "$cli_reads_fraction_mode" ]] && set_config_value reads_fraction_mode "$cli_reads_fraction_mode"
        [[ -n "$cli_reads_cutoff" ]] && set_config_value reads_cutoff "$cli_reads_cutoff"
        [[ -n "$cli_slope_cutoff" ]] && set_config_value slope_cutoff "$cli_slope_cutoff"
        ;;
    image)
        output_path=$(dirname "$input_abs")
        if $input_from_cli || [[ -z ${tissue_positions_file:-} || -z ${fullres_image_path:-} || -z ${frame_mask_file:-} ]]; then
            set_config_value image_path "$input_abs"
            set_config_value image_result_path "$output_path"
            set_config_value frame_mask_file "$output_path/mask.png"
            set_config_value tissue_positions_file "$output_path/tissue_positions.tsv.gz"
            set_config_value tissue_mask_file "$output_path/tissue_mask.png"
            set_config_value fullres_image_path "$input_abs"
        fi
        ;;
esac

if $chip_from_cli; then
    set_config_value chip "$selected_chip"
fi

source "$config_abs"

if [[ -n "$selected_chip" ]]; then
    apply_chip_preset "$selected_chip" || {
        echo "Error: unsupported chip '$selected_chip'." >&2
        echo "Valid chips: $(chip_preset_names_csv)." >&2
        exit 1
    }
fi
export chip x_spots_number y_spots_number length_spot interval
export barcode_a_whitelist_path barcode_b_whitelist_path

case "$step" in
    mrna)
        script="$STEP_SCRIPT_DIR/mrna.sh"
        cpus=$sbatch_mrna_cpus; partition=$sbatch_mrna_partition
        memory=$sbatch_mrna_mem; walltime=$sbatch_mrna_time
        ;;
    saturation)
        script="$STEP_SCRIPT_DIR/saturation.sh"
        cpus=$sbatch_mrna_cpus; partition=$sbatch_mrna_partition
        memory=$sbatch_mrna_mem; walltime=$sbatch_mrna_time
        ;;
    darlin)
        script="$STEP_SCRIPT_DIR/darlin.sh"
        cpus=$sbatch_darlin_cpus; partition=$sbatch_darlin_partition
        memory=$sbatch_darlin_mem; walltime=$sbatch_darlin_time
        ;;
    image)
        script="$STEP_SCRIPT_DIR/image.sh"
        cpus=$sbatch_image_cpus; partition=$sbatch_image_partition
        memory=$sbatch_image_mem; walltime=$sbatch_image_time
        ;;
esac

script_args=("$config_abs")
if [[ "$step" == saturation ]]; then
    script_args+=("$cli_saturation_fractions")
fi

if [[ ${execution_mode} == local ]]; then
    echo "Running $step locally"
    exec bash "$script" "${script_args[@]}"
fi

if [[ -z "$partition" ]]; then
    echo "Error: SLURM partition is required for step '$step'." >&2
    exit 1
fi

sbatch_args=(
    -J "${sbatch_job_name_prefix}_${step}"
    -c "$cpus"
    -p "$partition"
    --mem="$memory"
    --time="$walltime"
    -o "$sbatch_output"
    -e "$sbatch_error"
    --export="ALL,STEP_SCRIPT_DIR=$STEP_SCRIPT_DIR,REPO_DIR=$REPO_DIR,chip=$chip,x_spots_number=$x_spots_number,y_spots_number=$y_spots_number,length_spot=$length_spot,interval=$interval,barcode_a_whitelist_path=$barcode_a_whitelist_path,barcode_b_whitelist_path=$barcode_b_whitelist_path"
)
if [[ "${sbatch_requeue:-false}" =~ ^([Tt][Rr][Uu][Ee]|[Yy][Ee][Ss]|1)$ ]]; then
    sbatch_args+=(--requeue)
fi
echo "Submitting $step"
sbatch "${sbatch_args[@]}" "$script" "${script_args[@]}"
