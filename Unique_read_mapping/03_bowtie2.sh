#!/bin/bash
#SBATCH --array=0-81
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8000

# array currently set for data from publication
# e.g. #SBATCH --array=0-81 for 82 samples in array - 3x PBRM1, 3x SMARCA4, 3x CENPB, 3x H3K9me2, 3 x H3K9me3 and 8x IgG in parental RPE1,
# and 3x SMARCA4, 3x H3K9me2, 3x H3K9me3 and 6x IgG in PBRM1KO RPE1, 3x SMARCA4 in SMARCA4KO RPE1 - low and high salt for each


# Requires approx 64G memory and 24 hours running time (8 threads with 8G per core)

# Align reads to reference genome using Bowtie2 ##

# When performing CUT&RUN on human cells with yeast spike-in DNA, reads are first aligned to a combined human-yeast reference genome ##

cd ..

INPUT_FILES=(trimmed_reads/*_1_val_1.fq.gz)
base=$(basename -s _1_val_1.fq.gz ${INPUT_FILES[$SLURM_ARRAY_TASK_ID]})

bowtie2 -p 8 \
-x genomes/combined_t2tsc/CHM13-T2T_yeast/combined \
-1 trimmed_reads/${base}_1_val_1.fq.gz \
-2 trimmed_reads/${base}_2_val_2.fq.gz \
-S sam/${base}_t2tsc.sam &> sam/${base}_t2tsc.txt \
--local --very-sensitive-local --no-unal --no-mixed --no-discordant --dovetail --soft-clipped-unmapped-tlen \
--non-deterministic --phred33 -I 50 -X 1500
