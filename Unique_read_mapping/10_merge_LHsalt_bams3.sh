#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=4000

# array currently set for data from publication

# WT_IgG has 8 replicates so need to do last 2 separately

## If performing CUT&RUN with salt fractionation, the low and high salt libraries for each sample can be merged for downstream analysis ##

## Requires approx 1GB memory and 1 hour running time

# Only run using data analysed from publication - available on GEO

samtools merge bam/merged_bam/WT_IgG_LH_r7.t2t.NM4.sorted.final.bam \
 bam/WT_IgG_L_r7.t2t.NM4.sorted.final.bam bam/WT_IgG_H_r7.t2t.NM4.sorted.final.bam

samtools merge bam/merged_bam/WT_IgG_LH_r8.t2t.NM4.sorted.final.bam \
 bam/WT_IgG_L_r8.t2t.NM4.sorted.final.bam bam/WT_IgG_H_r8.t2t.NM4.sorted.final.bam
