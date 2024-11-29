#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=4000

# currently set for data from publication - only necessary if merging low and high salt bams where one has had duplicates removed

# Replace LH_salt merged bams with the low_salt_removed_duplicates bam files merged to high_salt bams for two specific samples

# Requires approx 1GB memory and 1 hour running time


rm -f bam/merged_bam/WT_IgG_LH_r3.t2t.NM4.sorted.final.bam

rm -f bam/merged_bam/WT_SMARCA4_LH_r1.t2t.NM4.sorted.final.bam

samtools merge bam/merged_bam/WT_IgG_LrmDupH_r3.t2t.NM4.sorted.final.bam \
 bam/WT_IgG_L_r3.t2t.NM4.sorted.final.rmDup.bam bam/WT_IgG_H_r3.t2t.NM4.sorted.final.bam

samtools merge bam/merged_bam/WT_SMARCA4_LrmDupH_r1.t2t.NM4.sorted.final.bam \
 bam/WT_SMARCA4_L_r1.t2t.NM4.sorted.final.rmDup.bam bam/WT_SMARCA4_H_r1.t2t.NM4.sorted.final.bam

