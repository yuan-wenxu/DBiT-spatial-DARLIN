#!/bin/bash

# Show help
show_help() {
    cat << EOF
Usage: $0 -f <fastq_path> [-o <output_dir>] [OPTIONS]

Process mRNA sequencing data with preprocessing, STAR alignment, and quality control.

Required Arguments:
  -f, --fastq_path <dir>            Directory containing R1/R2 fastq files

Output Options:
  -o, --output_path <dir>           Output directory path (optional; default: parent directory of fastq_path)

Preprocessing Options:
  -w, --whitelist <path>            Path to barcode whitelist file

  Advanced options:
    --compression_level <num>        Compression level for gzip (default: 6)
    --linker1 <seq>                  Linker 1 sequence (default: GTGGCCGATGTTTCGCATCGGCGTACGACT)
    --linker2 <seq>                  Linker 2 sequence (default: ATCCACGTGCTTGAGAGGCCAGAGCATTCG)
    --mm_rate <float>                Mismatch rate for linker sequences (default: 0.05)
    --bc_max_dist <num>              Maximum distance for barcode correction (default: 1)
    --scratch <path>                 Path to scratch directory for intermediate files (optional)

STAR Alignment Options:
  --genome_dir <path>               Path to STAR genome index
  --star_threads <num>              Number of threads for STAR (default: 10)
  --solo_cb_start <num>             Start position of cell barcode (default: 1)
  --solo_cb_len <num>               Length of cell barcode (default: 16)
  --solo_umi_start <num>            Start position of UMI (default: 17)
  --solo_umi_len <num>              Length of UMI (default: 10)

mRNA QC Options:
  --umi_min <num>                   Minimum UMI count per spot (default: 900)
  --gene_min <num>                  Minimum gene count per spot (default: 300)
  --min_cells <num>                 Minimum number of cells per gene (default: 3)
  --x_spots_number <num>            Number of spots in x direction (default: 50)
  --y_spots_number <num>            Number of spots in y direction (default: 50 )
  --length_spot <num>               Length of each spot in pixels (default: 20)
  --interval <num>                  Interval between spots in pixels (default: 20)
  --pixel_length <float>            Length of each pixel in microns (default: 0.294)

Pixi environment options:
  --pixi_env <name>                   Name of the Pixi environment to use (optional; default: dbit)

Other Options:
  -h, --help                        Show this help message and exit

Examples:
  # Basic usage
  $0 -f /path/to/fastq

  # Write results to a different directory
  $0 -f /path/to/fastq -o /path/to/output

EOF
}

# Set default values

# Preprocessing Advanced options
compression_level=${compression_level:-6}
linker1=${linker1:-GTGGCCGATGTTTCGCATCGGCGTACGACT}
linker2=${linker2:-ATCCACGTGCTTGAGAGGCCAGAGCATTCG}
mm_rate=${mm_rate:-0.05}
bc_max_dist=${bc_max_dist:-1}
scratch=${scratch:-}

# STAR Alignment Options
star_threads=${star_threads:-10}
soloCBstart=${soloCBstart:-1}
soloCBlen=${soloCBlen:-16}
soloUMIstart=${soloUMIstart:-17}
soloUMIlen=${soloUMIlen:-10}

# mRNA QC Options
umi_min=${umi_min:-900}
gene_min=${gene_min:-300}
min_cells=${min_cells:-3}
x_spots_number=${x_spots_number:-50}
y_spots_number=${y_spots_number:-50}
length_spot=${length_spot:-20}
interval=${interval:-20}
pixel_length=${pixel_length:-0.294}

# Pixi environment options
pixi_env=${pixi_env:-dbit}

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd) || exit 1
PYTHON_DIR="$SCRIPT_DIR/python"

normalize_dir_path() {
    local path="$1"
    while [[ "$path" != "/" && "$path" == */ ]]; do
        path="${path%/}"
    done
    printf '%s\n' "$path"
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

while getopts "f:o:w:h" opt; do
    case $opt in
        f) fastq_path=$OPTARG ;;
        o) output_path=$OPTARG ;;
        w) whitelist_path=$OPTARG ;;
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
        # Preprocessing Advanced options
        --compression_level) compression_level=$2; shift 2 ;;
        --linker1) linker1=$2; shift 2 ;;
        --linker2) linker2=$2; shift 2 ;;
        --mm_rate) mm_rate=$2; shift 2 ;;
        --bc_max_dist) bc_max_dist=$2; shift 2 ;;
        --scratch) scratch=$2; shift 2 ;;
        # STAR Alignment Options
        --genome_dir) genome_dir=$2; shift 2 ;;
        --star_threads) star_threads=$2; shift 2 ;;
        --solo_cb_start) soloCBstart=$2; shift 2 ;;
        --solo_cb_len) soloCBlen=$2; shift 2 ;;
        --solo_umi_start) soloUMIstart=$2; shift 2 ;;
        --solo_umi_len) soloUMIlen=$2; shift 2 ;;
        # mRNA QC Options
        --umi_min) umi_min=$2; shift 2 ;;
        --gene_min) gene_min=$2; shift 2 ;;
        --min_cells) min_cells=$2; shift 2 ;;
        --x_spots_number) x_spots_number=$2; shift 2 ;;
        --y_spots_number) y_spots_number=$2; shift 2 ;;
        --length_spot) length_spot=$2; shift 2 ;;
        --interval) interval=$2; shift 2 ;;
        --pixel_length) pixel_length=$2; shift 2 ;;
        --pixi_env) pixi_env=$2; shift 2 ;;
        # Other Options
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

fastq_path=$(normalize_dir_path "$fastq_path")
if [ -n "$output_path" ]; then
    output_path=$(normalize_dir_path "$output_path")
fi
if [ -n "$scratch" ]; then
    scratch=$(normalize_dir_path "$scratch")
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
    bam_file="$orig_output_path/results/$sample_name/Aligned.sortedByCoord.out.bam"

    pre_file=("$orig_output_path"/*barcode/${sample_name}_bc_match_R1.fq.gz)
    pre_done=false
    star_done=false
    if ((${#pre_file[@]} > 0)) && [ -f "$log_file" ]; then pre_done=true; fi
    if [ -f "$bam_file" ]; then star_done=true; fi

    use_scratch=false
    if [ -n "$scratch" ]; then
        use_scratch=true
        scratch_sample="$scratch/$sample_name"
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

        pixi run -e "$pixi_env" python "$PYTHON_DIR/preprocess.py" \
            -r1 "$step1_r1" -r2 "$step1_r2" \
            -o "$step1_out" -s "$sample_name" \
            -b1 "$whitelist_path" -b2 "$whitelist_path" \
            -cl "$compression_level" \
            -l1 "$linker1" -l2 "$linker2" -m "$mm_rate" \
            -bc_max_dist "$bc_max_dist" &> "$step1_log"

        # Keep step1 outputs in original output path for future skip checks.
        if $use_scratch; then
            cp "$step1_log" "$log_file"
            cp -r "$scratch_output"/*barcode "$orig_output_path"/ 2>/dev/null || true
        fi
    fi

    # Resolve step2 input from existing/preprocessed files in original output path.
    pre_file=("$orig_output_path"/*barcode/${sample_name}_bc_match_R1.fq.gz)
    if ((${#pre_file[@]} > 0)); then
        pre_r1="${pre_file[0]}"
        pre_r2="${pre_r1%_R1.fq.gz}_R2.fq.gz"
        tmp_path="$(dirname "$pre_r1")"
    else
        echo "Error: missing preprocess outputs for $sample_name, skip STAR."
        if $use_scratch; then rm -rf "$scratch_sample"; fi
        continue
    fi

    # Step 2: STAR (skip if already done)
    if $star_done; then
        echo "Step2 STAR already done for $sample_name, skipping..."
    else
        if $use_scratch; then
            mkdir -p "$scratch_input" "$scratch_output"
            cp "$pre_r1" "$pre_r2" "$scratch_input/"
            star_input="$scratch_input"
            star_results="$scratch_output/results/$sample_name"
        else
            star_input="$tmp_path"
            star_results="$orig_output_path/results/$sample_name"
        fi

        mkdir -p "$star_results"
        pixi run -e "$pixi_env" STAR \
            --runMode alignReads \
            --runThreadN "$star_threads" \
            --genomeDir "$genome_dir" \
            --readFilesCommand zcat \
            --readFilesIn "$star_input/${sample_name}_bc_match_R2.fq.gz" "$star_input/${sample_name}_bc_match_R1.fq.gz" \
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
            --soloCellFilter EmptyDrops_CR \
            --soloFeatures GeneFull \
            --bamRemoveDuplicatesType UniqueIdentical \
            --quantMode GeneCounts &> "$star_results/STAR.log"

        # Step3 is local only, so copy step2 results back first when using scratch.
        if $use_scratch; then
            mkdir -p "$orig_output_path/results/$sample_name"
            cp -r "$star_results"/* "$orig_output_path/results/$sample_name"/
        fi
    fi

    # Step 3: always run locally, no skip and no scratch.
    final_results="$orig_output_path/results/$sample_name"
    mkdir -p "$final_results"
    pixi run -e "$pixi_env" python "$PYTHON_DIR/mrna.py" -f "$final_results/Solo.out" -w "$whitelist_path" \
        -umi_min "$umi_min" -gene_min "$gene_min" -min_cells "$min_cells" \
        --x_spots_number "$x_spots_number" --y_spots_number "$y_spots_number" \
        --length_spot "$length_spot" --interval "$interval" \
        --pixel_length "$pixel_length" &> "$final_results/Solo.out/qc.log"

    if $use_scratch; then
        rm -rf "$scratch_sample"
    fi
done
