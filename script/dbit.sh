#!/bin/bash
set -o pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd) || exit 1
QC_REPO_DIR=$(cd "$SCRIPT_DIR/.." && pwd) || exit 1
QC_SCRIPT_DIR="$SCRIPT_DIR/Quality_Control"
export QC_SCRIPT_DIR QC_REPO_DIR

show_help() {
    cat <<EOF
Usage: $0 <step> [options]

Steps:
  mrna          Process transcriptome FASTQs and run spatial mRNA QC
  amplicon      Process DARLIN amplicon FASTQs
  image         Segment a registered image and count cells
  plot          Final step: generate cell-filtered plots after image

Run '$0 <step> -h' to show parameters for one step.
EOF
}

show_mrna_help() {
    cat <<EOF
Usage: $0 mrna --input <fastq_dir> --config <file> [options]

Required:
  --input <path>    Transcriptome FASTQ directory
  --config <file>   Per-dataset configuration file

Optional:
  --chip <name>     50-50, 50-20, or 100-20; required only before stored
  --umi-min <int>   Minimum UMI count per spot (default: 900)
  --gene-min <int>  Minimum gene count per spot (default: 300)
  --min-cell <int>  Minimum cells per gene (default: 3)
EOF
}

show_amplicon_help() {
    cat <<EOF
Usage: $0 amplicon --input <fastq_dir> --config <file> [options]

Required:
  --input <path>    Amplicon FASTQ directory
  --config <file>   Per-dataset configuration file

Optional:
  --chip <name>                                   50-50, 50-20, or 100-20; required only before stored
  --initial-reads-cutoff <int>                    Reads cutoff for initial filtering (default: 100)
  --major-fraction-threshold-molecule <float>     Reads fraction threshold for major molecule (default: 0.8)
  --reads-fraction-mode <sum|max>                 Mode for calculating reads fraction (default: sum)
  --reads-cutoff <int>                            Reads cutoff for final filtering (default: 10)
  --slope-cutoff <float>                          Slope cutoff for final filtering (default: 10)
EOF
}

show_image_help() {
    cat <<EOF
Usage: $0 image --input <image> --config <file> [options]

Required:
  --input <path>          Registered input image
  --config <file>         Per-dataset configuration file
  --orientation <mode>    normal, horizontal, vertical, or rotate
  --swap-xy <bool>        True or False

Optional:
  --chip <name>           50-50, 50-20, or 100-20; required only before stored
EOF
}

show_plot_help() {
    cat <<EOF
Usage: $0 plot --config <file> [--chip <name>]

Required:
  --config <file>   Configuration populated by earlier data steps

Optional:
  --chip <name>     50-50, 50-20, or 100-20; required only before stored
EOF
}

show_step_help() {
    case "$1" in
        mrna) show_mrna_help ;;
        amplicon) show_amplicon_help ;;
        image) show_image_help ;;
        plot) show_plot_help ;;
    esac
}

if [[ $# -eq 0 || ${1:-} == -h || ${1:-} == --help ]]; then show_help; exit 0; fi

step=$1
if [[ ${2:-} == -h || ${2:-} == --help ]]; then
    case "$step" in
        mrna|amplicon|image|plot) show_step_help "$step"; exit 0 ;;
    esac
fi
shift

case "$step" in
    mrna|amplicon|image|plot) ;;
    *)
        echo "Error: unsupported step '$step'." >&2
        echo "Valid steps: mrna, amplicon, image, plot." >&2
        exit 1
        ;;
esac

input_path=""
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
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input) input_path=$2; shift 2 ;;
        --config) config_file=$2; shift 2 ;;
        --chip) selected_chip=$2; chip_from_cli=true; shift 2 ;;
        --umi-min) cli_umi_min=$2; shift 2 ;;
        --gene-min) cli_gene_min=$2; shift 2 ;;
        --min-cell) cli_min_cell=$2; shift 2 ;;
        --initial-reads-cutoff) cli_initial_reads_cutoff=$2; shift 2 ;;
        --major-fraction-threshold-molecule) cli_major_fraction_threshold_molecule=$2; shift 2 ;;
        --reads-fraction-mode) cli_reads_fraction_mode=$2; shift 2 ;;
        --reads-cutoff) cli_reads_cutoff=$2; shift 2 ;;
        --slope-cutoff) cli_slope_cutoff=$2; shift 2 ;;
        --orientation) cli_orientation=$2; shift 2 ;;
        --swap-xy) cli_swap_xy=$2; shift 2 ;;
        *) shift ;;
    esac
done

if [[ "$step" != mrna && ( -n "$cli_umi_min" || -n "$cli_gene_min" || -n "$cli_min_cell" ) ]]; then
    echo "Error: --umi-min, --gene-min, and --min-cell can only be used with the mrna step." >&2
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
if [[ "$step" == image && ( -z "$cli_orientation" || -z "$cli_swap_xy" ) ]]; then
    echo "Error: image requires both --orientation and --swap-xy." >&2
    exit 1
fi
if [[ "$step" == plot && -n "$input_path" ]]; then
    echo "Error: --input cannot be used with the plot step; paths are read from the config." >&2
    exit 1
fi

if [[ ! -f "$config_file" ]]; then
    echo "Error: config file not found: $config_file" >&2
    exit 1
fi
config_abs=$(realpath "$config_file")

set_config_value() {
    local key=$1
    local value=$2
    printf '%s=%q\n' "$key" "$value" >> "$config_abs"
}

if ! $chip_from_cli; then
    chip=""
    source "$config_abs"
    selected_chip=${chip:-}
fi
if [[ -z "$selected_chip" ]]; then
    echo "Error: --chip is required the first time; no chip is stored in $config_abs." >&2
    exit 1
fi
case "$selected_chip" in
    50-50|50-20|100-20) ;;
    *)
        echo "Error: unsupported chip '$selected_chip'." >&2
        echo "Valid chips: 50-50, 50-20, 100-20." >&2
        exit 1
        ;;
esac

if [[ "$step" != plot ]]; then
    input_abs=$(realpath -m "$input_path")
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
        set_config_value mrna_fastq_path "$input_abs"
        set_config_value mrna_output_path "$output_path"
        set_config_value mrna_dir "$output_path/results/$sample_name/Solo.out/GeneFull"
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
        set_config_value amplicon_fastq_path "$input_abs"
        set_config_value amplicon_output_path "$output_path"
        set_config_value amp_dir "$output_path/results/$nonlocus_sample_name"
        [[ -n "$cli_initial_reads_cutoff" ]] && set_config_value initial_reads_cutoff "$cli_initial_reads_cutoff"
        [[ -n "$cli_major_fraction_threshold_molecule" ]] && set_config_value major_fraction_threshold_molecule "$cli_major_fraction_threshold_molecule"
        [[ -n "$cli_reads_fraction_mode" ]] && set_config_value reads_fraction_mode "$cli_reads_fraction_mode"
        [[ -n "$cli_reads_cutoff" ]] && set_config_value reads_cutoff "$cli_reads_cutoff"
        [[ -n "$cli_slope_cutoff" ]] && set_config_value slope_cutoff "$cli_slope_cutoff"
        ;;
    image)
        output_path=$(dirname "$input_abs")
        set_config_value image_path "$input_abs"
        set_config_value image_result_path "$output_path"
        set_config_value cell_number_file "$output_path/filtered_results.csv"
        set_config_value gray_path "$output_path/gray.png"
        [[ -n "$cli_orientation" ]] && set_config_value orientation "$cli_orientation"
        [[ -n "$cli_swap_xy" ]] && set_config_value swap_xy "$cli_swap_xy"
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
        whitelist_path="$QC_REPO_DIR/docs/barcodes/barcodes.tsv"
        ;;
    50-20)
        chip=50-20; x_spots_number=50; y_spots_number=50
        length_spot=20; interval=20
        whitelist_path="$QC_REPO_DIR/docs/barcodes/barcodes.tsv"
        ;;
    100-20)
        chip=100-20; x_spots_number=100; y_spots_number=100
        length_spot=20; interval=20
        whitelist_path="$QC_REPO_DIR/docs/barcodes/barcodes100.tsv"
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
esac

if [[ ${execution_mode} == local ]]; then
    echo "Running $step locally"
    exec bash "$script" "$config_abs"
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
    --export="ALL,QC_SCRIPT_DIR=$QC_SCRIPT_DIR,QC_REPO_DIR=$QC_REPO_DIR,chip=$chip,x_spots_number=$x_spots_number,y_spots_number=$y_spots_number,length_spot=$length_spot,interval=$interval,whitelist_path=$whitelist_path"
)
if [[ "${sbatch_requeue:-false}" =~ ^([Tt][Rr][Uu][Ee]|[Yy][Ee][Ss]|1)$ ]]; then
    sbatch_args+=(--requeue)
fi
echo "Submitting $step"
sbatch "${sbatch_args[@]}" "$script" "$config_abs"
