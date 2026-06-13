#!/bin/bash

# Show help message
show_help() {
    cat << EOF
Usage: $0 -f <fastq_path> [-o <output_path>] [OPTIONS]

Process amplification sequencing data with preprocessing, DARLIN correction, and visualization.

Required Arguments:
  -f, --fastq_path <dir>            Directory containing R1/R2 fastq files

Output Options:
  -o, --output_path <dir>           Output directory path (optional; default: parent directory of fastq_path)

Preprocessing Options:
  -w, --whitelist <path>            Path to barcode whitelist file
  -c, --cutadapt <bool>             Perform cutadapt trimming (True/False) (default: True)

  Advanced options:
    --cores <num>                     Number of cores for cutadapt and barcode extraction (default: 8)
    --preprocess_batch_size <num>     Read pairs per barcode extraction worker batch (default: 50000)
    --base_quality <num>              Base quality score threshold (default: 10)
    --compression_level <num>         Compression level for gzip barcode FASTQ output (default: 6)
    --gzip_output <bool>              Gzip barcode FASTQ output (default: false). false uses more disk space but can speed up preprocessing.
    --gzip_after_preprocess <bool>    If gzip_output is false, compress barcode FASTQs after preprocessing with pigz/gzip (default: true)
    --linker1 <seq>                   Linker 1 sequence (default: GTGGCCGATGTTTCGCATCGGCGTACGACT)
    --linker2 <seq>                   Linker 2 sequence (default: ATCCACGTGCTTGAGAGGCCAGAGCATTCG)
    --mm_rate <float>                 Mismatch rate for linker sequences (default: 0.05)
    --scratch <path>                  Path to scratch directory for intermediate files (optional)

DARLIN Correction Options:
  --darlin <bool>                    Whether lineage barcode sequences are available (default: True)
  --sb_len <num>                     Length of concatenated spot barcode (default: 16)
  --ub_len <num>                     Length of UMI barcode (default: 10)
  --umi_hd_threshold <num>           Hamming-distance threshold for UMI correction within each SR (default: 1)
  --min_lb_len <num>                 Minimum lineage barcode length (default: 20)
  --initial_reads_cutoff <num>       Minimum reads per raw LB/SB/UB molecule (default: 100)
  --lb_error_rate <float>            Lineage barcode correction error rate (default: 0.01)
  --lb_min_hd <num>                  Minimum HD threshold for lineage barcode correction (default: 1)
  --major_fraction_threshold <float> Minimum major LR fraction per SR/UR (default: 0.8)
  --reads_fraction_mode <sum|max>    Denominator for major LR filtering (default: sum)
  --reads_cutoff <num>               Minimum reads per final SR/UR/LR row (default: 10)
  --slope_cutoff <num>               Minimum reads/UMIs per SR (default: 10)

Pixi environment options:
  --pixi_env <name>                   Name of the Pixi environment to use (optional; default: dbit)
  --pixi_env_dir <path>               Directory containing pixi.toml (optional; default: repository root)

Other Options:
  -h, --help                        Show this help message and exit

Examples:
  # Basic usage
  $0 -f /path/to/fastq -c True

  # Write results to a different directory
  $0 -f /path/to/fastq -o /path/to/output -c True

EOF
}

# Set default values
# Preprocessing Options
cutadapt=${cutadapt:-True}

# Preprocessing Advanced options
cores=${cores:-8}
preprocess_batch_size=${preprocess_batch_size:-50000}
base_quality=${base_quality:-10}
compression_level=${compression_level:-6}
gzip_output=${gzip_output:-false}
gzip_after_preprocess=${gzip_after_preprocess:-true}
linker1=${linker1:-GTGGCCGATGTTTCGCATCGGCGTACGACT}
linker2=${linker2:-ATCCACGTGCTTGAGAGGCCAGAGCATTCG}
mm_rate=${mm_rate:-0.05}
scratch=${scratch:-}

# DARLIN Correction Options
darlin=${darlin:-True}
sb_len=${sb_len:-16}
ub_len=${ub_len:-10}
umi_hd_threshold=${umi_hd_threshold:-1}
min_lb_len=${min_lb_len:-20}
initial_reads_cutoff=${initial_reads_cutoff:-100}
lb_error_rate=${lb_error_rate:-0.01}
lb_min_hd=${lb_min_hd:-1}
major_fraction_threshold_molecule=${major_fraction_threshold_molecule:-0.8}
reads_fraction_mode=${reads_fraction_mode:-sum}
reads_cutoff=${reads_cutoff:-10}
slope_cutoff=${slope_cutoff:-10}

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd) || exit 1
PYTHON_DIR="$SCRIPT_DIR/python"
whitelist_path="$SCRIPT_DIR/../../docs/barcodes/barcodes.tsv"

# Pixi environment options
pixi_env=${pixi_env:-dbit}
pixi_env_dir=${pixi_env_dir:-$(cd "$SCRIPT_DIR/../.." && pwd)}

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

# Short options
short_args=()
while [[ $# -gt 0 ]]; do
    case $1 in
        --*) break ;;
        *) short_args+=("$1"); shift ;;
    esac
done
long_args=("$@")

set -- "${short_args[@]}"
OPTIND=1

while getopts "f:o:w:c:h" opt; do
    case $opt in
        f) fastq_path=$OPTARG ;;
        o) output_path=$OPTARG ;;
        w) whitelist_path=$OPTARG ;;
        c) cutadapt=$OPTARG  ;;
        h) show_help; exit 0 ;;
        ?) echo "Invalid option: -$OPTARG" >&2
            echo "Use -h or --help for usage information" >&2
            exit 1 ;;
    esac
done

# Long options
set -- "${long_args[@]}"
while [[ $# -gt 0 ]]; do
    case $1 in
        # Required Arguments
        --fastq_path) fastq_path=$2; shift 2 ;;
        # Output Options
        --output_path) output_path=$2; shift 2 ;;
        # Preprocessing Options
        --whitelist) whitelist_path=$2; shift 2 ;;
        --cutadapt) cutadapt=$2; shift 2 ;;
        # Preprocessing Advanced options
        --cores) cores=$2; shift 2 ;;
        --preprocess_batch_size) preprocess_batch_size=$2; shift 2 ;;
        --base_quality) base_quality=$2; shift 2 ;;
        --compression_level) compression_level=$2; shift 2 ;;
        --gzip_output) gzip_output=$2; shift 2 ;;
        --gzip_after_preprocess) gzip_after_preprocess=$2; shift 2 ;;
        --linker1) linker1=$2; shift 2 ;;
        --linker2) linker2=$2; shift 2 ;;
        --mm_rate) mm_rate=$2; shift 2 ;;
        --scratch) scratch=$2; shift 2 ;;
        # DARLIN Correction Options
        --darlin) darlin=$2; shift 2 ;;
        --sb_len) sb_len=$2; shift 2 ;;
        --ub_len) ub_len=$2; shift 2 ;;
        --umi_hd_threshold) umi_hd_threshold=$2; shift 2 ;;
        --min_lb_len) min_lb_len=$2; shift 2 ;;
        --initial_reads_cutoff) initial_reads_cutoff=$2; shift 2 ;;
        --lb_error_rate) lb_error_rate=$2; shift 2 ;;
        --lb_min_hd) lb_min_hd=$2; shift 2 ;;
        --major_fraction_threshold_molecule) major_fraction_threshold_molecule=$2; shift 2 ;;
        --reads_fraction_mode) reads_fraction_mode=$2; shift 2 ;;
        --reads_cutoff) reads_cutoff=$2; shift 2 ;;
        --slope_cutoff) slope_cutoff=$2; shift 2 ;;
        --pixi_env) pixi_env=$2; shift 2 ;;
        --pixi_env_dir) pixi_env_dir=$2; shift 2 ;;
        --help) show_help; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

# Check required arguments
if [ -z "$fastq_path" ]; then
    echo "Error: -f (fastq_path) is required" >&2
    echo "Use -h or --help for usage information" >&2
    exit 1
fi

if [ -z "$whitelist_path" ]; then
    echo "Error: -w (whitelist) is required" >&2
    echo "Use -h or --help for usage information" >&2
    exit 1
fi

if [[ "$reads_fraction_mode" != "sum" && "$reads_fraction_mode" != "max" ]]; then
    echo "Error: --reads_fraction_mode must be 'sum' or 'max'" >&2
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
    scratch_input="$scratch/amp/input"
    scratch_output="$scratch/amp/output"
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
    #locus=${sample_name: -2}
    r2=$file_path/$sample_name"_R2.fq.gz"
    nonlocus_sample_name=$(echo "$sample_name" | sed 's/\(-CA\|-RA\|-TA\)//')

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
        -cb "false" &> "$output_path/${sample_name}_preprocess.log"

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
        -o "$results" -d "$darlin" \
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
        --slope-cutoff "$slope_cutoff" &> "$results/dbit.log"

    run_pixi python "$PYTHON_DIR/plot/heatmap.py" \
        -f "$results/final.csv" \
        -w "$whitelist_path" \
        -o "$results"
done

if [ -n "$scratch" ]; then
    cp -r "$scratch_output"/* "$orig_output_path"/
    rm -rf "$scratch/amp"
fi
