#!/bin/bash
# DBiT-spatial-DARLIN QC pipeline configuration

# Basic configuration
genome_dir=/path/to/STAR_genome_index
scratch=                         # Optional temporary directory
execution_mode=hpc               # hpc or local
pixel_length=0.294
seed=42

# mRNA
mrna_cores=10
compression_level=6
gzip_output=false
gzip_after_preprocess=true
linker1=GTGGCCGATGTTTCGCATCGGCGTACGACT
linker2=ATCCACGTGCTTGAGAGGCCAGAGCATTCG
mm_rate=0.05
preprocess_batch_size=50000
star_threads=${mrna_cores}
soloCBstart=1
soloCBlen=16
soloUMIstart=17
soloUMIlen=10

# Saturation analysis
saturation_fractions=0.01,0.02,0.05,0.1,0.2,0.5

# DARLIN
darlin_cores=10
base_quality=10
cutadapt=true
sb_len=${soloCBlen}
ub_len=${soloUMIlen}
umi_hd_threshold=1
min_lb_len=20
lb_error_rate=0.01
lb_min_hd=1

# SLURM resources (used when execution_mode=hpc)
sbatch_job_name_prefix=dbit
sbatch_output=%x_%j.out
sbatch_error=%x_%j.err
sbatch_requeue=true

sbatch_mrna_cpus=${mrna_cores}
sbatch_mrna_partition=
sbatch_mrna_mem=64G
sbatch_mrna_time=24:00:00

sbatch_darlin_cpus=${darlin_cores}
sbatch_darlin_partition=
sbatch_darlin_mem=64G
sbatch_darlin_time=24:00:00

sbatch_image_cpus=8
sbatch_image_partition=
sbatch_image_mem=32G
sbatch_image_time=12:00:00

# File paths and chip
