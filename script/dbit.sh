#!/bin/bash
set -o pipefail

SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || exit 1
SCRIPT_DIR=$(cd "$(dirname "$SCRIPT_PATH")" && pwd) || exit 1
REPO_DIR=$(cd "$SCRIPT_DIR/.." && pwd) || exit 1
QC_SCRIPT_DIR="$SCRIPT_DIR/Quality_Control"
LR_SCRIPT_DIR="$SCRIPT_DIR/Clone_Analysis"
PROGRAM_NAME=$(basename "$0")
export QC_SCRIPT_DIR REPO_DIR LR_SCRIPT_DIR

show_help() {
    cat <<EOF
Usage: $PROGRAM_NAME <step> [options]

Steps:
  mrna          Process transcriptome FASTQs and run spatial mRNA QC
  amplicon      Process DARLIN amplicon FASTQs
  image         Segment a registered image and count cells
  plot          Generate tissue-filtered plots after image
  clone         Filter and plot clone-analysis results

Run '$PROGRAM_NAME <step> -h' to show parameters for one step.
EOF
}

show_mrna_help() {
    cat <<EOF
Usage: $PROGRAM_NAME mrna --config <file> [options]

Required:
  --config <file>   Per-dataset configuration file

Optional:
  --input <path>         Transcriptome FASTQ directory; required only before stored
  --chip <name>          50-50, 50-20, or 100-20; required only before stored
  --genome-dir <path>    STAR genome index directory; overrides the config value
  --umi-min <int>        Non-negative minimum UMI count per spot (default: 900)
  --gene-min <int>       Non-negative minimum gene count per spot (default: 300)
  --min-cell <int>       Positive minimum cells per gene (default: 3)
EOF
}

show_amplicon_help() {
    cat <<EOF
Usage: $PROGRAM_NAME amplicon --config <file> [options]

Required:
  --config <file>   Per-dataset configuration file

Optional:
  --input <path>                                  Amplicon FASTQ directory; required only before stored
  --chip <name>                                   50-50, 50-20, or 100-20; required only before stored
  --initial-reads-cutoff <int>                    Non-negative reads cutoff for initial filtering (default: 100)
  --major-fraction-threshold-molecule <float>     Major-molecule reads fraction from 0 to 1 (default: 0.8)
  --reads-fraction-mode <sum|max>                 Mode for calculating reads fraction (default: sum)
  --reads-cutoff <int>                            Non-negative reads cutoff for final filtering (default: 10)
  --slope-cutoff <float>                          Non-negative slope cutoff for final filtering (default: 10)
EOF
}

show_image_help() {
    cat <<EOF
Usage: $PROGRAM_NAME image --config <file> [options]

Required:
  --config <file>         Per-dataset configuration file

Optional:
  --input <path>          Registered input image; required only before stored
  --orientation <mode>    normal, horizontal, vertical, or rotate; required only before stored
  --swap-xy <bool>        True or False (case-insensitive); required only before stored
  --chip <name>           50-50, 50-20, or 100-20; required only before stored

90-degree rotation combinations:
  --orientation horizontal --swap-xy True    90 degrees counterclockwise
  --orientation vertical   --swap-xy True    90 degrees clockwise
EOF
}

show_plot_help() {
    cat <<EOF
Usage: $PROGRAM_NAME plot --config <file> [--chip <name>]

Required:
  --config <file>   Configuration populated by earlier data steps

Optional:
  --chip <name>     50-50, 50-20, or 100-20; required only before stored
EOF
}

show_clone_help() {
    cat <<EOF
Usage: $PROGRAM_NAME clone --config <file> [options]

Required:
  --config <file>     Configuration populated by earlier data steps

Optional:
  --labels <list>     Comma-separated labels (default: CA,RA,TA)
  --top-n <int>       Positive number of LR plots per label (default: 10)
EOF
}

show_step_help() {
    case "$1" in
        mrna) show_mrna_help ;;
        amplicon) show_amplicon_help ;;
        image) show_image_help ;;
        plot) show_plot_help ;;
        clone) show_clone_help ;;
    esac
}

if [[ $# -eq 0 || ${1:-} == -h || ${1:-} == --help ]]; then show_help; exit 0; fi

step=$1
if [[ ${2:-} == -h || ${2:-} == --help ]]; then
    case "$step" in
        mrna|amplicon|image|plot|clone) show_step_help "$step"; exit 0 ;;
    esac
fi
shift

case "$step" in
    mrna|amplicon|image|plot|clone) ;;
    *)
        echo "Error: unsupported step '$step'." >&2
        echo "Valid steps: mrna, amplicon, image, plot, clone." >&2
        exit 1
        ;;
esac

input_path=""
input_from_cli=false
config_file=""
selected_chip=""
chip_from_cli=false
cli_umi_min=""
cli_gene_min=""
cli_min_cell=""
cli_genome_dir=""
cli_initial_reads_cutoff=""
cli_major_fraction_threshold_molecule=""
cli_reads_fraction_mode=""
cli_reads_cutoff=""
cli_slope_cutoff=""
cli_orientation=""
cli_swap_xy=""
cli_clone_labels=""
cli_top_n=""

require_option_value() {
    if [[ $# -lt 2 || $2 == --* ]]; then
        echo "Error: option '$1' requires a value." >&2
        exit 1
    fi
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input) require_option_value "$@"; input_path=$2; input_from_cli=true; shift 2 ;;
        --config) require_option_value "$@"; config_file=$2; shift 2 ;;
        --chip) require_option_value "$@"; selected_chip=$2; chip_from_cli=true; shift 2 ;;
        --genome-dir) require_option_value "$@"; cli_genome_dir=$2; shift 2 ;;
        --umi-min) require_option_value "$@"; cli_umi_min=$2; shift 2 ;;
        --gene-min) require_option_value "$@"; cli_gene_min=$2; shift 2 ;;
        --min-cell) require_option_value "$@"; cli_min_cell=$2; shift 2 ;;
        --initial-reads-cutoff) require_option_value "$@"; cli_initial_reads_cutoff=$2; shift 2 ;;
        --major-fraction-threshold-molecule) require_option_value "$@"; cli_major_fraction_threshold_molecule=$2; shift 2 ;;
        --reads-fraction-mode) require_option_value "$@"; cli_reads_fraction_mode=$2; shift 2 ;;
        --reads-cutoff) require_option_value "$@"; cli_reads_cutoff=$2; shift 2 ;;
        --slope-cutoff) require_option_value "$@"; cli_slope_cutoff=$2; shift 2 ;;
        --orientation) require_option_value "$@"; cli_orientation=$2; shift 2 ;;
        --swap-xy) require_option_value "$@"; cli_swap_xy=$2; shift 2 ;;
        --labels) require_option_value "$@"; cli_clone_labels=$2; shift 2 ;;
        --top-n) require_option_value "$@"; cli_top_n=$2; shift 2 ;;
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

if [[ "$step" != mrna && ( -n "$cli_umi_min" || -n "$cli_gene_min" || -n "$cli_min_cell" ) ]]; then
    echo "Error: --umi-min, --gene-min, and --min-cell can only be used with the mrna step." >&2
    exit 1
fi
if [[ "$step" != mrna && -n "$cli_genome_dir" ]]; then
    echo "Error: --genome-dir can only be used with the mrna step." >&2
    exit 1
fi
if [[ "$step" != amplicon && ( -n "$cli_initial_reads_cutoff" || -n "$cli_major_fraction_threshold_molecule" || -n "$cli_reads_fraction_mode" || -n "$cli_reads_cutoff" || -n "$cli_slope_cutoff" ) ]]; then
    echo "Error: amplicon filtering options can only be used with the amplicon step." >&2
    exit 1
fi
if [[ "$step" != image && ( -n "$cli_orientation" || -n "$cli_swap_xy" ) ]]; then
    echo "Error: --orientation and --swap-xy can only be used with the image step." >&2
    exit 1
fi
if [[ "$step" != clone && ( -n "$cli_clone_labels" || -n "$cli_top_n" ) ]]; then
    echo "Error: --labels and --top-n can only be used with the clone step." >&2
    exit 1
fi
if [[ "$step" == plot && -n "$input_path" ]]; then
    echo "Error: --input cannot be used with the plot step; paths are read from the config." >&2
    exit 1
fi
if [[ "$step" == clone && -n "$input_path" ]]; then
    echo "Error: --input cannot be used with the clone step; amp_dir is read from the config." >&2
    exit 1
fi

validate_nonnegative_integer --umi-min "$cli_umi_min"
validate_nonnegative_integer --gene-min "$cli_gene_min"
validate_positive_integer --min-cell "$cli_min_cell"
validate_nonnegative_integer --initial-reads-cutoff "$cli_initial_reads_cutoff"
validate_fraction --major-fraction-threshold-molecule "$cli_major_fraction_threshold_molecule"
validate_nonnegative_integer --reads-cutoff "$cli_reads_cutoff"
validate_nonnegative_number --slope-cutoff "$cli_slope_cutoff"
validate_positive_integer --top-n "$cli_top_n"
if [[ -n "$cli_reads_fraction_mode" ]]; then
    case "$cli_reads_fraction_mode" in
        sum|max) ;;
        *)
            echo "Error: --reads-fraction-mode must be sum or max; got '$cli_reads_fraction_mode'." >&2
            exit 1
            ;;
    esac
fi

if [[ ! -f "$config_file" ]]; then
    echo "Error: config file not found: $config_file" >&2
    exit 1
fi
config_abs=$(realpath "$config_file")

# Load stored values first. Command-line values below take precedence and are
# appended to the config so that later runs can reuse them.
source "$config_abs"

set_config_value() {
    local key=$1
    local value=$2
    printf '%s=%q\n' "$key" "$value" >> "$config_abs"
}

if ! $chip_from_cli; then selected_chip=${chip:-}; fi
if [[ "$step" != clone && -z "$selected_chip" ]]; then
    echo "Error: --chip is required the first time; no chip is stored in $config_abs." >&2
    exit 1
fi
case "$selected_chip" in
    50-50|50-20|100-20|"") ;;
    *)
        echo "Error: unsupported chip '$selected_chip'." >&2
        echo "Valid chips: 50-50, 50-20, 100-20." >&2
        exit 1
        ;;
esac

if [[ "$step" != plot && -z "$input_path" ]]; then
    case "$step" in
        mrna) input_path=${mrna_fastq_path:-} ;;
        amplicon) input_path=${amplicon_fastq_path:-} ;;
        image) input_path=${image_path:-} ;;
        clone) input_path=${amp_dir} ;;
    esac
fi
if [[ "$step" != plot && "$step" != clone && -z "$input_path" ]]; then
    echo "Error: --input is required the first time; no input path for '$step' is stored in $config_abs." >&2
    exit 1
fi
if [[ "$step" == mrna ]]; then
    effective_genome_dir=${cli_genome_dir:-${genome_dir:-}}
    if [[ -n "$cli_genome_dir" ]]; then
        effective_genome_dir=$(realpath -m "$cli_genome_dir")
    fi
    if [[ -z "$effective_genome_dir" ]]; then
        echo "Error: --genome-dir is required the first time; no genome_dir is stored in $config_abs." >&2
        exit 1
    fi
    if [[ ! -d "$effective_genome_dir" ]]; then
        echo "Error: genome directory does not exist: $effective_genome_dir" >&2
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
if [[ "$step" != plot ]]; then
    input_abs=$(realpath -m "$input_path")
fi
if [[ "$step" == image && ! -f "$input_abs" ]]; then
    echo "Error: image file does not exist: $input_abs" >&2
    exit 1
fi
if [[ "$step" == clone && ! -d "$input_abs" ]]; then
    echo "Error: clone input directory does not exist: $input_abs" >&2
    exit 1
fi

case "$step" in
    mrna)
        first_r1=$(find "$input_abs" -maxdepth 1 -type f -name '*_R1.fq.gz' -print -quit)
        if [[ -z "$first_r1" ]]; then
            echo "Error: no *_R1.fq.gz file found in $input_abs" >&2
            exit 1
        fi
        sample_name=$(basename "$first_r1" _R1.fq.gz)
        output_path=$(dirname "$input_abs")
        if $input_from_cli; then
            set_config_value mrna_fastq_path "$input_abs"
            set_config_value mrna_output_path "$output_path"
            set_config_value mrna_dir "$output_path/results/$sample_name/Solo.out/GeneFull"
            set_config_value cluster_csv "$output_path/results/$sample_name/Solo.out/GeneFull/raw/data_tissuefiltered.csv"
        fi
        [[ -n "$cli_genome_dir" ]] && set_config_value genome_dir "$effective_genome_dir"
        [[ -n "$cli_umi_min" ]] && set_config_value umi_min "$cli_umi_min"
        [[ -n "$cli_gene_min" ]] && set_config_value gene_min "$cli_gene_min"
        [[ -n "$cli_min_cell" ]] && set_config_value min_cells "$cli_min_cell"
        ;;
    amplicon)
        first_r1=$(find "$input_abs" -maxdepth 1 -type f -name '*_R1.fq.gz' -print -quit)
        if [[ -z "$first_r1" ]]; then
            echo "Error: no *_R1.fq.gz file found in $input_abs" >&2
            exit 1
        fi
        sample_name=$(basename "$first_r1" _R1.fq.gz)
        nonlocus_sample_name=$(echo "$sample_name" | sed 's/\(-CA\|-RA\|-TA\|_CA\|_RA\|_TA\)//')
        output_path=$(dirname "$input_abs")
        if $input_from_cli; then
            set_config_value amplicon_fastq_path "$input_abs"
            set_config_value amplicon_output_path "$output_path"
            set_config_value amp_dir "$output_path/results/$nonlocus_sample_name"
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
    clone)
        [[ -n "$cli_clone_labels" ]] && set_config_value clone_labels "$cli_clone_labels"
        [[ -n "$cli_top_n" ]] && set_config_value clone_top_n "$cli_top_n"
        ;;
esac

if $chip_from_cli; then
    set_config_value chip "$selected_chip"
fi

source "$config_abs"

# Chip presets are intentionally defined only in this submission script.
case "$selected_chip" in
    50-50)
        chip=50-50; x_spots_number=50; y_spots_number=50
        length_spot=50; interval=50
        whitelist_path="$REPO_DIR/docs/barcodes/barcodes.tsv"
        ;;
    50-20)
        chip=50-20; x_spots_number=50; y_spots_number=50
        length_spot=20; interval=20
        whitelist_path="$REPO_DIR/docs/barcodes/barcodes.tsv"
        ;;
    100-20)
        chip=100-20; x_spots_number=100; y_spots_number=100
        length_spot=20; interval=20
        whitelist_path="$REPO_DIR/docs/barcodes/barcodes100.tsv"
        ;;
esac
export chip x_spots_number y_spots_number length_spot interval whitelist_path

case "$step" in
    mrna)
        script="$QC_SCRIPT_DIR/mrna.sh"
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
    plot)
        script="$QC_SCRIPT_DIR/plot.sh"
        cpus=$sbatch_plot_cpus; partition=$sbatch_plot_partition
        memory=$sbatch_plot_mem; walltime=$sbatch_plot_time
        ;;
    clone)
        script="$LR_SCRIPT_DIR/top_lr_pipeline.sh"
        cpus=${sbatch_clone_cpus}; partition=${sbatch_clone_partition}
        memory=${sbatch_clone_mem}; walltime=${sbatch_clone_time}
        ;;
esac

script_args=("$config_abs")

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
    --export="ALL,QC_SCRIPT_DIR=$QC_SCRIPT_DIR,LR_SCRIPT_DIR=$LR_SCRIPT_DIR,REPO_DIR=$REPO_DIR,chip=$chip,x_spots_number=$x_spots_number,y_spots_number=$y_spots_number,length_spot=$length_spot,interval=$interval,whitelist_path=$whitelist_path"
)
if [[ "${sbatch_requeue:-false}" =~ ^([Tt][Rr][Uu][Ee]|[Yy][Ee][Ss]|1)$ ]]; then
    sbatch_args+=(--requeue)
fi
echo "Submitting $step"
sbatch "${sbatch_args[@]}" "$script" "${script_args[@]}"
