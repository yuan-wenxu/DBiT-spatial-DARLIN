#!/bin/bash

# Show help message
show_help() {
    cat << EOF
Usage: $0 -f <reads_dir> -o <output_dir> -d <cutadapt> [OPTIONS]

Process amplification sequencing data with preprocessing, DARLIN correction, and visualization.

Required Arguments:
  -f, --reads_dir <dir>             Directory name containing R1/R2 fastq files (relative to output_path) (default: fastq)
  -o, --output_path <dir>           Output directory path

Preprocessing Options:
  -w, --whitelist <path>            Path to barcode whitelist file
  -c, --cutadapt <bool>             Perform cutadapt trimming (True/False) (default: False)

  Advanced options:
    --cores <num>                    Number of cores for cutadapt  (default: 8)
    --base_quality <num>             Base quality score threshold (default: 10)
    --compression_level <num>        Compression level for gzip (default: 1)
    --linker1 <seq>                  Linker 1 sequence (default: CAAGCGTTGGCTTCTCGCATCT)
    --linker2 <seq>                  Linker 2 sequence (default: ATCCACGTGCTTGAGAGGCCAGAGCATTCG)
    --tn5 <seq>                      Tn5 sequence (default: GTGGCCGATGTTTCGCATCGGCGTACGACT)
    --mm_rate <float>                Mismatch rate for linker sequences (default: 0.05)
    --use_linker1 <bool>             Use linker 1 for barcode correction (True/False) (default: False)
    --bc_max_dist <num>              Maximum distance for barcode correction (default: 1)
    --scratch <path>                 Path to scratch directory for intermediate files (optional)

DARLIN Correction Options:
  --umi_hd_threshold <num>            Edit-distance threshold for UMI clustering within each SR (default: 1)
  --lb_error_rate <float>             Per-base error rate used to derive LB edit-distance threshold (default: 0.02)
  --major_fraction_threshold <float>  Per (CR, UR) group: keep LR with reads fraction >= this value (default: 0.8)
  --reads_cutoff <num>                Only keep (SR, UR, LR) with supported reads >= this value (default: 10)
  --slope_cutoff <num>                Only keep SR (spots) with k = reads/UMIs >= this value (default: 10)

Other Options:
  -h, --help                        Show this help message and exit

Examples:
  # Basic usage
  $0 -f reads -o /path/to/output -c True

EOF
}

# Set default values
reads_dir=${reads_dir:-fastq}
# Preprocessing Options
cutadapt=${cutadapt:-False}

# Preprocessing Advanced options
cores=${cores:-8}
base_quality=${base_quality:-10}
compression_level=${compression_level:-1}
linker1=${linker1:-CAAGCGTTGGCTTCTCGCATCT}
linker2=${linker2:-ATCCACGTGCTTGAGAGGCCAGAGCATTCG}
tn5=${tn5:-GTGGCCGATGTTTCGCATCGGCGTACGACT}
mm_rate=${mm_rate:-0.05}
use_linker1=${use_linker1:-False}
bc_max_dist=${bc_max_dist:-1}
scratch=${scratch:-}

# DARLIN Correction Options
umi_hd_threshold=${umi_hd_threshold:-1}
lb_error_rate=${lb_error_rate:-0.02}
major_fraction_threshold_molecule=${major_fraction_threshold_molecule:-0.8}
reads_cutoff=${reads_cutoff:-10}
slope_cutoff=${slope_cutoff:-10}

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
        f) reads_dir=$OPTARG ;;
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
        # Preprocessing Advanced options
        --cores) cores=$2; shift 2 ;;
        --base_quality) base_quality=$2; shift 2 ;;
        --compression_level) compression_level=$2; shift 2 ;;
        --linker1) linker1=$2; shift 2 ;;
        --linker2) linker2=$2; shift 2 ;;
        --tn5) tn5=$2; shift 2 ;;
        --mm_rate) mm_rate=$2; shift 2 ;;
        --use_linker1) use_linker1=$2; shift 2 ;;
        --bc_max_dist) bc_max_dist=$2; shift 2 ;;
        --scratch) scratch=$2; shift 2 ;;
        # DARLIN Correction Options
        --umi_hd_threshold) umi_hd_threshold=$2; shift 2 ;;
        --lb_error_rate) lb_error_rate=$2; shift 2 ;;
        --major_fraction_threshold) major_fraction_threshold_molecule=$2; shift 2 ;;
        --reads_cutoff) reads_cutoff=$2; shift 2 ;;
        --slope_cutoff) slope_cutoff=$2; shift 2 ;;
        --help) show_help; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

# Check required arguments
if [ -z "$reads_dir" ] || [ -z "$output_path" ]; then
    echo "Error: -f (reads_dir) and -o (output_path) are required" >&2
    echo "Use -h or --help for usage information" >&2
    exit 1
fi

if [ -n "$scratch" ]; then
    scratch_input="$scratch/input"
    scratch_output="$scratch/output"
    mkdir -p "$scratch_input" "$scratch_output"
    cp -r "$output_path/$reads_dir"/* "$scratch_input/"
    orig_output_path="$output_path"
    file_path="$scratch_input"
    output_path="$scratch_output"
else
    file_path="$output_path/$reads_dir"
fi

for r1 in $file_path/*_R1.fq.gz; do
    sample_name=$(basename $r1 | sed 's/_R1.fq.gz//')
    locus=$(echo "$sample_name" | grep -oE 'CA|RA|TA')
    #locus=${sample_name: -2}
    r2=$file_path/$sample_name"_R2.fq.gz"
    nonlocus_sample_name=$(echo "$sample_name" | sed 's/\(-CA\|-RA\|-TA\)//')

    # Cutadapt and extract UMI and barcode
    python ./python/preprocess.py \
        -r1 "$r1" -r2 "$r2" \
        -o "$output_path" -s "$sample_name" \
        -b1 "$whitelist_path" -b2 "$whitelist_path" \
        -l "$locus" -c "$cores" -q "$base_quality" \
        -cl "$compression_level" -cut "$cutadapt" \
        -l1 "$linker1" -l2 "$linker2" -tn5 "$tn5" \
        -m "$mm_rate" -ul1 "$use_linker1" \
        -bc_max_dist "$bc_max_dist" &> "$output_path/${sample_name}_preprocess.log"

    tmp_path=$(tail -n 1 $output_path/${sample_name}_preprocess.log)

    results=$output_path/results/$nonlocus_sample_name/$locus
    mkdir -p "$results"

    python ./python/amplicon.py \
        -bu "$tmp_path/${sample_name}_bc_match_R1.fq.gz" \
        -dr "$tmp_path/${sample_name}_bc_match_R2.fq.gz" \
        -o "$results" -d "$cutadapt" \
        --umi_hd_threshold "$umi_hd_threshold" \
        --lb_error_rate "$lb_error_rate" \
        --major_fraction_threshold_molecule "$major_fraction_threshold_molecule" \
        --reads_cutoff "$reads_cutoff" \
        --slope_cutoff "$slope_cutoff" &> "$results/dbit.log"

    python ./python/plot/heatmap.py \
        -f "$results/final.csv" \
        -w "$whitelist_path" \
        -o "$results"
done

if [ -n "$scratch" ]; then
    cp -r "$scratch_output"/* "$orig_output_path"/
    rm -rf "$scratch"
fi