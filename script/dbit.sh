#!/bin/bash
set -o pipefail

SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || exit 1
SCRIPT_DIR=$(cd "$(dirname "$SCRIPT_PATH")" && pwd) || exit 1
REPO_DIR=$(cd "$SCRIPT_DIR/.." && pwd) || exit 1
START_DIR=$(pwd -P) || exit 1
QC_SCRIPT_DIR="$SCRIPT_DIR/Quality_Control"
LR_SCRIPT_DIR="$SCRIPT_DIR/Clone_Analysis"
DOMAIN_SCRIPT_DIR="$SCRIPT_DIR/Domain_Analysis"
SATURATION_SCRIPT_DIR="$SCRIPT_DIR/Saturation"
PROGRAM_NAME=$(basename "$0")
CHIP_FILE="$REPO_DIR/config/chip.sh"
export QC_SCRIPT_DIR REPO_DIR LR_SCRIPT_DIR DOMAIN_SCRIPT_DIR SATURATION_SCRIPT_DIR

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
  amplicon      Process DARLIN amplicon FASTQs
  image         Segment a registered image and count cells
  filter        Apply the tissue mask and generate filtered plots
  domain        Run RCTD deconvolution and BANKSY domain clustering
  clone         Filter and plot clone-analysis results

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

show_amplicon_help() {
    cat <<EOF
Usage: $PROGRAM_NAME amplicon [--config <file>] [options]

Optional:
  --config <file>                                 Configuration file (default: ./dbit.config.sh)
  --input <path>                                  Amplicon FASTQ directory; required only before stored
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
  --input <path>          Registered input image; required only before stored
  --orientation <mode>    normal, horizontal, vertical, or rotate; required only before stored
  --swap-xy <bool>        True or False (case-insensitive); required only before stored

90-degree rotation combinations:
  --orientation horizontal --swap-xy True    90 degrees counterclockwise
  --orientation vertical   --swap-xy True    90 degrees clockwise
EOF
}

show_filter_help() {
    cat <<EOF
Usage: $PROGRAM_NAME filter [--config <file>]

Optional:
  --config <file>   Configuration file (default: ./dbit.config.sh)
EOF
}

show_clone_help() {
    cat <<EOF
Usage: $PROGRAM_NAME clone [--config <file>] [options]

Optional:
  --config <file>      Configuration file (default: ./dbit.config.sh)
  --labels <list>      Comma-separated labels (default: CA,RA,TA)
  --top-n <int>        Positive number of LR plots per label (default: 10)
  --rotate <degrees>   Clockwise grid rotation for display: 0, 90, 180, or 270
EOF
}

show_domain_help() {
    cat <<EOF
Usage: $PROGRAM_NAME domain [--config <file>] [options]

Optional:
  --config <file>        Configuration file (default: ./dbit.config.sh)
  --rotate <degrees>     Clockwise grid rotation for display: 0, 90, 180, or 270
EOF
}

show_step_help() {
    case "$1" in
        init) show_init_help ;;
        mrna) show_mrna_help ;;
        saturation) show_saturation_help ;;
        amplicon) show_amplicon_help ;;
        image) show_image_help ;;
        filter) show_filter_help ;;
        domain) show_domain_help ;;
        clone) show_clone_help ;;
    esac
}

if [[ $# -eq 0 || ${1:-} == -h || ${1:-} == --help ]]; then show_help; exit 0; fi

step=$1
if [[ ${2:-} == -h || ${2:-} == --help ]]; then
    case "$step" in
        init|mrna|saturation|amplicon|image|filter|domain|clone) show_step_help "$step"; exit 0 ;;
    esac
fi
shift

case "$step" in
    init|mrna|saturation|amplicon|image|filter|domain|clone) ;;
    *)
        echo "Error: unsupported step '$step'." >&2
        echo "Valid steps: init, mrna, saturation, amplicon, image, filter, domain, clone." >&2
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
cli_orientation=""
cli_swap_xy=""
cli_clone_labels=""
cli_top_n=""
cli_rotate=""
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
        --config) require_step_option "$1" mrna saturation amplicon image filter domain clone; require_option_value "$@"; config_file=$2; shift 2 ;;
        --input) require_step_option "$1" mrna amplicon image; require_option_value "$@"; input_path=$2; input_from_cli=true; shift 2 ;;
        --chip) require_step_option "$1" mrna amplicon; require_option_value "$@"; selected_chip=$2; chip_from_cli=true; shift 2 ;;
        --umi-min) require_step_option "$1" mrna; require_option_value "$@"; cli_umi_min=$2; shift 2 ;;
        --gene-min) require_step_option "$1" mrna; require_option_value "$@"; cli_gene_min=$2; shift 2 ;;
        --min-cell) require_step_option "$1" mrna; require_option_value "$@"; cli_min_cell=$2; shift 2 ;;
        --fractions) require_step_option "$1" saturation; require_option_value "$@"; cli_saturation_fractions=$2; shift 2 ;;
        --initial-reads-cutoff) require_step_option "$1" amplicon; require_option_value "$@"; cli_initial_reads_cutoff=$2; shift 2 ;;
        --major-fraction-threshold-molecule) require_step_option "$1" amplicon; require_option_value "$@"; cli_major_fraction_threshold_molecule=$2; shift 2 ;;
        --reads-fraction-mode) require_step_option "$1" amplicon; require_option_value "$@"; cli_reads_fraction_mode=$2; shift 2 ;;
        --reads-cutoff) require_step_option "$1" amplicon; require_option_value "$@"; cli_reads_cutoff=$2; shift 2 ;;
        --slope-cutoff) require_step_option "$1" amplicon; require_option_value "$@"; cli_slope_cutoff=$2; shift 2 ;;
        --orientation) require_step_option "$1" image; require_option_value "$@"; cli_orientation=$2; shift 2 ;;
        --swap-xy) require_step_option "$1" image; require_option_value "$@"; cli_swap_xy=$2; shift 2 ;;
        --labels) require_step_option "$1" clone; require_option_value "$@"; cli_clone_labels=$2; shift 2 ;;
        --top-n) require_step_option "$1" clone; require_option_value "$@"; cli_top_n=$2; shift 2 ;;
        --rotate) require_step_option "$1" domain clone; require_option_value "$@"; cli_rotate=$2; shift 2 ;;
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

validate_nonnegative_integer --umi-min "$cli_umi_min"
validate_nonnegative_integer --gene-min "$cli_gene_min"
validate_positive_integer --min-cell "$cli_min_cell"
validate_nonnegative_integer --initial-reads-cutoff "$cli_initial_reads_cutoff"
validate_fraction --major-fraction-threshold-molecule "$cli_major_fraction_threshold_molecule"
validate_nonnegative_integer --reads-cutoff "$cli_reads_cutoff"
validate_nonnegative_number --slope-cutoff "$cli_slope_cutoff"
validate_positive_integer --top-n "$cli_top_n"
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
if [[ -n "$cli_rotate" ]]; then
    case "$cli_rotate" in
        0|90|180|270) ;;
        *)
            echo "Error: --rotate must be 0, 90, 180, or 270; got '$cli_rotate'." >&2
            exit 1
            ;;
    esac
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
    amplicon) stored_input=${amplicon_fastq_path:-} ;;
    image) stored_input=${image_path:-} ;;
esac
[[ -n "$input_path" ]] || input_path=$stored_input
if [[ "$step" =~ ^(mrna|saturation|amplicon|image)$ && -z "$input_path" ]]; then
    if [[ "$step" == saturation ]]; then
        echo "Error: mrna_fastq_path is not stored in $config_abs; run 'dbit mrna --input <fastq_dir>' first." >&2
    else
        echo "Error: --input is required the first time; no input path for '$step' is stored in $config_abs." >&2
    fi
    exit 1
fi
if [[ "$step" == mrna ]]; then
    if [[ -z ${genome_dir:-} ]]; then
        echo "Error: genome_dir must be set in $config_abs." >&2
        exit 1
    fi
fi
if [[ "$step" == image ]]; then
    effective_orientation=${cli_orientation:-${orientation:-}}
    effective_swap_xy=${cli_swap_xy:-${swap_xy:-}}
    if [[ -z "$effective_orientation" || -z "$effective_swap_xy" ]]; then
        echo "Error: --orientation and --swap-xy are required the first time; no values are stored in $config_abs." >&2
        exit 1
    fi
    case "$effective_orientation" in
        normal|horizontal|vertical|rotate) ;;
        *)
            echo "Error: --orientation must be normal, horizontal, vertical, or rotate; got '$effective_orientation'." >&2
            exit 1
            ;;
    esac
    case "${effective_swap_xy,,}" in
        true) effective_swap_xy=True ;;
        false) effective_swap_xy=False ;;
        *)
            echo "Error: --swap-xy must be True or False; got '$effective_swap_xy'." >&2
            exit 1
            ;;
    esac
fi
if [[ "$step" == mrna || "$step" == saturation || "$step" == amplicon || "$step" == image ]]; then
    input_abs=$(realpath -m "$input_path")
fi

case "$step" in
    mrna)
        first_r1=$(find "$input_abs" -maxdepth 1 -type f -name '*_R1.fq.gz' -print -quit)
        if [[ -z "$first_r1" ]]; then
            echo "Error: no *_R1.fq.gz file found in $input_abs" >&2
            exit 1
        fi
        output_path=$(dirname "$input_abs")
        if $input_from_cli; then
            set_config_value mrna_fastq_path "$input_abs"
            set_config_value mrna_output_path "$output_path"
            set_config_value mrna_dir "$output_path/results/Solo.out/GeneFull"
            set_config_value cluster_csv "$output_path/results/Solo.out/GeneFull/raw/data_tissuefiltered.csv"
        fi
        [[ -n "$cli_umi_min" ]] && set_config_value umi_min "$cli_umi_min"
        [[ -n "$cli_gene_min" ]] && set_config_value gene_min "$cli_gene_min"
        [[ -n "$cli_min_cell" ]] && set_config_value min_cells "$cli_min_cell"
        ;;
    saturation)
        first_r1=$(find "$input_abs" -maxdepth 1 -type f -name '*_R1.fq.gz' -print -quit)
        if [[ -z "$first_r1" ]]; then
            echo "Error: no *_R1.fq.gz file found in $input_abs" >&2
            exit 1
        fi
        ;;
    amplicon)
        first_r1=$(find "$input_abs" -maxdepth 1 -type f -name '*_R1.fq.gz' -print -quit)
        if [[ -z "$first_r1" ]]; then
            echo "Error: no *_R1.fq.gz file found in $input_abs" >&2
            exit 1
        fi
        output_path=$(dirname "$input_abs")
        if $input_from_cli; then
            set_config_value amplicon_fastq_path "$input_abs"
            set_config_value amplicon_output_path "$output_path"
            set_config_value amp_dir "$output_path/results"
        fi
        [[ -n "$cli_initial_reads_cutoff" ]] && set_config_value initial_reads_cutoff "$cli_initial_reads_cutoff"
        [[ -n "$cli_major_fraction_threshold_molecule" ]] && set_config_value major_fraction_threshold_molecule "$cli_major_fraction_threshold_molecule"
        [[ -n "$cli_reads_fraction_mode" ]] && set_config_value reads_fraction_mode "$cli_reads_fraction_mode"
        [[ -n "$cli_reads_cutoff" ]] && set_config_value reads_cutoff "$cli_reads_cutoff"
        [[ -n "$cli_slope_cutoff" ]] && set_config_value slope_cutoff "$cli_slope_cutoff"
        ;;
    image)
        output_path=$(dirname "$input_abs")
        if $input_from_cli; then
            set_config_value image_path "$input_abs"
            set_config_value image_result_path "$output_path"
            set_config_value cell_number_file "$output_path/filtered_results.csv"
            set_config_value tissue_mask_file "$output_path/tissue_mask.png"
            set_config_value gray_path "$output_path/gray.png"
        fi
        [[ -n "$cli_orientation" ]] && set_config_value orientation "$cli_orientation"
        [[ -n "$cli_swap_xy" ]] && set_config_value swap_xy "$effective_swap_xy"
        ;;
    domain)
        [[ -n "$cli_rotate" ]] && set_config_value rotate "$cli_rotate"
        ;;
    clone)
        [[ -n "$cli_clone_labels" ]] && set_config_value clone_labels "$cli_clone_labels"
        [[ -n "$cli_top_n" ]] && set_config_value clone_top_n "$cli_top_n"
        [[ -n "$cli_rotate" ]] && set_config_value rotate "$cli_rotate"
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
export chip x_spots_number y_spots_number length_spot interval whitelist_path

case "$step" in
    mrna)
        script="$QC_SCRIPT_DIR/mrna.sh"
        cpus=$sbatch_mrna_cpus; partition=$sbatch_mrna_partition
        memory=$sbatch_mrna_mem; walltime=$sbatch_mrna_time
        ;;
    saturation)
        script="$SATURATION_SCRIPT_DIR/saturation.sh"
        cpus=$sbatch_mrna_cpus; partition=$sbatch_mrna_partition
        memory=$sbatch_mrna_mem; walltime=$sbatch_mrna_time
        ;;
    amplicon)
        script="$QC_SCRIPT_DIR/amplicon.sh"
        cpus=$sbatch_amplicon_cpus; partition=$sbatch_amplicon_partition
        memory=$sbatch_amplicon_mem; walltime=$sbatch_amplicon_time
        ;;
    image)
        script="$QC_SCRIPT_DIR/image.sh"
        cpus=$sbatch_image_cpus; partition=$sbatch_image_partition
        memory=$sbatch_image_mem; walltime=$sbatch_image_time
        ;;
    filter)
        script="$QC_SCRIPT_DIR/filter.sh"
        cpus=$sbatch_filter_cpus; partition=$sbatch_filter_partition
        memory=$sbatch_filter_mem; walltime=$sbatch_filter_time
        ;;
    domain)
        script="$DOMAIN_SCRIPT_DIR/domain.sh"
        cpus=${sbatch_domain_cpus}; partition=${sbatch_domain_partition}
        memory=${sbatch_domain_mem}; walltime=${sbatch_domain_time}
        ;;
    clone)
        script="$LR_SCRIPT_DIR/clone.sh"
        cpus=${sbatch_clone_cpus}; partition=${sbatch_clone_partition}
        memory=${sbatch_clone_mem}; walltime=${sbatch_clone_time}
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
    --export="ALL,QC_SCRIPT_DIR=$QC_SCRIPT_DIR,LR_SCRIPT_DIR=$LR_SCRIPT_DIR,DOMAIN_SCRIPT_DIR=$DOMAIN_SCRIPT_DIR,SATURATION_SCRIPT_DIR=$SATURATION_SCRIPT_DIR,REPO_DIR=$REPO_DIR,chip=$chip,x_spots_number=$x_spots_number,y_spots_number=$y_spots_number,length_spot=$length_spot,interval=$interval,whitelist_path=$whitelist_path"
)
if [[ "${sbatch_requeue:-false}" =~ ^([Tt][Rr][Uu][Ee]|[Yy][Ee][Ss]|1)$ ]]; then
    sbatch_args+=(--requeue)
fi
echo "Submitting $step"
sbatch "${sbatch_args[@]}" "$script" "${script_args[@]}"
