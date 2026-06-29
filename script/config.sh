#!/bin/bash
# DBiT-spatial-DARLIN QC pipeline configuration
# Run with: bash dbit.sh mrna /path/to/config.sh --chip 50-20

# =============================================================================
# Paths (edit these first)
# =============================================================================
mrna_fastq_path=/path/to/sample/transcriptome/fastq
mrna_output_path=/path/to/sample/transcriptome
amplicon_fastq_path=/path/to/sample/amplicon/fastq
amplicon_output_path=/path/to/sample/amplicon
image_path=/path/to/sample/image/align.png
image_result_path=/path/to/sample/image

# Inputs used by the final cell-filtered plotting step. These can point to
# outputs produced by the three steps above.
cell_number_file=/path/to/sample/image/filtered_results.csv
mrna_dir=/path/to/sample/transcriptome/results/sample/Solo.out/GeneFull
amp_dir=/path/to/sample/amplicon/results/sample
gray_path=/path/to/sample/image/gray.png

genome_dir=/path/to/STAR_genome_index
scratch=                         # Optional temporary directory

# =============================================================================
# Execution mode
# =============================================================================
execution_mode=hpc               # hpc or local

# =============================================================================
# Shared spatial parameter (chip presets are defined only in dbit.sh)
# =============================================================================
pixel_length=0.294
orientation=normal               # normal, horizontal, vertical, or rotate
swap_xy=False

# =============================================================================
# mRNA
# =============================================================================
mrna_cores=10
compression_level=6
gzip_output=false
gzip_after_preprocess=true
linker1=GTGGCCGATGTTTCGCATCGGCGTACGACT
linker2=ATCCACGTGCTTGAGAGGCCAGAGCATTCG
mm_rate=0.05
bc_max_dist=1
preprocess_batch_size=50000
star_threads=${mrna_cores}
soloCBstart=1
soloCBlen=16
soloUMIstart=17
soloUMIlen=10
umi_min=900
gene_min=300
min_cells=3

# =============================================================================
# Amplicon
# =============================================================================
amp_cores=10
base_quality=10
cutadapt=true
sb_len=${soloCBlen}
ub_len=${soloUMIlen}
umi_hd_threshold=1
min_lb_len=20
initial_reads_cutoff=100
lb_error_rate=0.01
lb_min_hd=1
major_fraction_threshold_molecule=0.8
reads_fraction_mode=sum           # sum or max
reads_cutoff=10
slope_cutoff=10

# =============================================================================
# Image
# =============================================================================
put_text=True
font_size=1
top_value=50
number_of_top_values=1500
model_name=2D_versatile_fluo
prob_thresh=0.5
nms_thresh=0.6
cutoff=100

# =============================================================================
# SLURM resources (used when execution_mode=hpc)
# =============================================================================
sbatch_job_name_prefix=dbit
sbatch_output=%x_%j.out
sbatch_error=%x_%j.err
sbatch_requeue=true

sbatch_mrna_cpus=${mrna_cores}
sbatch_mrna_partition=
sbatch_mrna_mem=64G
sbatch_mrna_time=24:00:00

sbatch_amplicon_cpus=${amp_cores}
sbatch_amplicon_partition=
sbatch_amplicon_mem=64G
sbatch_amplicon_time=24:00:00

sbatch_image_cpus=8
sbatch_image_partition=
sbatch_image_mem=32G
sbatch_image_time=12:00:00

sbatch_plot_cpus=1
sbatch_plot_partition=
sbatch_plot_mem=32G
sbatch_plot_time=04:00:00
