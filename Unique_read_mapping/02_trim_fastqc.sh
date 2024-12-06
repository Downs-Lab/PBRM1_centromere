#!/bin/bash
#SBATCH --array=0-81
#SBATCH --cpus-per-task=8
#SBATCH --time=6:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=4000

# array currently set for data from publication
# e.g. #SBATCH --array=0-81 for 82 samples in array - 3x PBRM1, 3x SMARCA4, 3x CENPB, 3x H3K9me2, 3 x H3K9me3 and 8x IgG in parental RPE1,
# and 3x SMARCA4, 3x H3K9me2, 3x H3K9me3 and 6x IgG in PBRM1KO RPE1, 3x SMARCA4 in SMARCA4KO RPE1 - low and high salt for each

# Requires approx 4G memory and <6 hours running time

cd ..

INPUT_FILES=(fastq/*_1.fq.gz)

base=$(basename -s _1.fq.gz ${INPUT_FILES[$SLURM_ARRAY_TASK_ID]})

# FASTQC quality control

fastqc fastq/${base}_1.fq.gz fastq/${base}_2.fq.gz -o fastqc

# trim files and get quality check output after (FASTQC)

trim_galore --trim-n --fastqc_args "--outdir fastqc/trim" --output_dir trimmed_reads --paired \
fastq/${base}_1.fq.gz \
fastq/${base}_2.fq.gz

