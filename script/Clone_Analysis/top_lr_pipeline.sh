#!/bin/bash

show_help() {
    cat << EOF
Usage: $0 <config_file>

Run the top-LR plotting pipeline.

Arguments:
  config_file   Per-dataset configuration file

Examples:
  $0 config.sh
EOF
}

SCRIPT_DIR=${CLONE_SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)} || exit 1
REPO_DIR=${REPO_DIR:-$(cd "$SCRIPT_DIR/../.." && pwd)} || exit 1
PYTHON_DIR="$SCRIPT_DIR/python"

if [[ ${1:-} == -h || ${1:-} == --help ]]; then show_help; exit 0; fi
if [[ $# -ne 1 ]]; then show_help >&2; exit 1; fi
config_file=$1
if [[ ! -f "$config_file" ]]; then
    echo "Error: config file not found: $config_file" >&2
    exit 1
fi

source "$config_file"
pixi_env=${pixi_env:-default}
pixi_env_dir=${pixi_env_dir:-$REPO_DIR}
input_dir=${amp_dir}
bank_dir=${bank_dir}
cluster_csv=${cluster_csv}
orientation=${orientation}
swap_xy=${swap_xy}
min_sequence_length=${min_sequence_length}
IFS=',' read -r -a labels <<< "${clone_labels:-CA,RA,TA}"
top_n=${clone_top_n:-10}


# Validate required config variables
if [[ -z "$input_dir" ]]; then
    echo "Error: amp_dir must be set in the config." >&2
    exit 1
fi
if [[ -z "$bank_dir" ]]; then
    echo "Error: bank_dir must be set in the config." >&2
    exit 1
fi
if [[ -z "$cluster_csv" ]]; then
    echo "Error: cluster_csv must be set in the config." >&2
    exit 1
fi
case "$orientation" in
    normal|horizontal|vertical|rotate) ;;
    *)
        echo "Error: orientation must be one of normal, horizontal, vertical, rotate; got '$orientation'." >&2
        exit 1
        ;;
esac
case "${swap_xy,,}" in
    true) swap_xy=True ;;
    false) swap_xy=False ;;
    *)
        echo "Error: swap_xy must be True or False; got '$swap_xy'." >&2
        exit 1
        ;;
esac


run_pixi() {
    (
        cd "$pixi_env_dir" || exit 1
        pixi run -e "$pixi_env" "$@"
    )
}

run_pixi python "$PYTHON_DIR/allele_bank_filter.py" \
    --input-dir "$input_dir" \
    --bank-dir "$bank_dir" \
    --min-sequence-length "$min_sequence_length" \
    --labels "${labels[@]}" || exit 1

run_pixi python "$PYTHON_DIR/top_lr_plot.py" \
    --input-dir "$input_dir" \
    --cluster-csv "$cluster_csv" \
    --top-n "$top_n" \
    --orientation "$orientation" \
    --labels "${labels[@]}" \
    $([ "$swap_xy" = True ] && printf %s --swap_xy) || exit 1
