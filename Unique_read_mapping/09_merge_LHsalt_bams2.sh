#!/bin/bash
#SBATCH --array=0-1
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=4000

# array currently set for data from publication

# WT_IgG and PBRM1KO_IgG have more than 3 replicates - so additional script with array for 2

## If performing CUT&RUN with salt fractionation, the low and high salt libraries for each sample can be merged for downstream analysis ##

## Requires approx 1GB memory and 2 hours running time

INPUT_FILES=(bam/*_L_r4.t2t.NM4.sorted.final.bam)

base=$(basename -s _L_r4.t2t.NM4.sorted.final.bam ${INPUT_FILES[$SLURM_ARRAY_TASK_ID]})


samtools merge bam/merged_bam/${base}_LH_r4.t2t.NM4.sorted.final.bam \
 bam/${base}_L_r4.t2t.NM4.sorted.final.bam bam/${base}_H_r4.t2t.NM4.sorted.final.bam

samtools merge bam/merged_bam/${base}_LH_r5.t2t.NM4.sorted.final.bam \
 bam/${base}_L_r5.t2t.NM4.sorted.final.bam bam/${base}_H_r5.t2t.NM4.sorted.final.bam

samtools merge bam/merged_bam/${base}_LH_r6.t2t.NM4.sorted.final.bam \
 bam/${base}_L_r6.t2t.NM4.sorted.final.bam bam/${base}_H_r6.t2t.NM4.sorted.final.bam



