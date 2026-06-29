#!/bin/bash
set -o pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd) || exit 1
QC_REPO_DIR=$(cd "$SCRIPT_DIR/.." && pwd) || exit 1
QC_SCRIPT_DIR="$SCRIPT_DIR/Quality_Control"
export QC_SCRIPT_DIR QC_REPO_DIR

show_help() {
    cat <<EOF
Usage: $0 <step> <config_file> --chip <chip_name>

Run or submit exactly one QC step.

Arguments:
  step          One of: mrna, amplicon, image, plot
  config_file   QC configuration file
  --chip        One of: 50-50, 50-20, 100-20

Example:
  $0 mrna config.sh --chip 50-20
EOF
}

if [[ ${1:-} == -h || ${1:-} == --help ]]; then show_help; exit 0; fi

step=$1
config_file=$2
selected_chip=$4

case "$step" in
    mrna|amplicon|image|plot) ;;
    *)
        echo "Error: unsupported step '$step'." >&2
        echo "Valid steps: mrna, amplicon, image, plot." >&2
        exit 1
        ;;
esac

# shellcheck disable=SC1090
source "$config_file"
config_abs=$(cd "$(dirname "$config_file")" && pwd)/$(basename "$config_file")

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
        script="$QC_SCRIPT_DIR/plot_cell_filtered.sh"
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
