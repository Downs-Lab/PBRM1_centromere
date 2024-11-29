#!/bin/bash
#SBATCH --array=0-10
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=4000

# array currently set for data from publication
# e.g. #SBATCH --array=0-10 for 11 different samples in array - PBRM1, SMARCA4, CENPB, H3K9me2, H3K9me3 and IgG in parental RPE1, and SMARCA4, H3K9me2, H3K9me3 and IgG in PBRM1KO RPE1,
# and SMARCA4 in SMARCA4KO RPE1 - replicates analysed together (up to 8 for WT_IgG)

## If performing CUT&RUN with salt fractionation, the low and high salt libraries for each sample can be merged for downstream analysis ##
## Requires approx 1GB memory and 2 hours running time

INPUT_FILES=(bam/*_L_r1.t2t.NM4.sorted.final.bam)

base=$(basename -s _L_r1.t2t.NM4.sorted.final.bam ${INPUT_FILES[$SLURM_ARRAY_TASK_ID]})


samtools merge bam/merged_bam/${base}_LH_r1.t2t.NM4.sorted.final.bam \
 bam/${base}_L_r1.t2t.NM4.sorted.final.bam bam/${base}_H_r1.t2t.NM4.sorted.final.bam

samtools merge bam/merged_bam/${base}_LH_r2.t2t.NM4.sorted.final.bam \
 bam/${base}_L_r2.t2t.NM4.sorted.final.bam bam/${base}_H_r2.t2t.NM4.sorted.final.bam

samtools merge bam/merged_bam/${base}_LH_r3.t2t.NM4.sorted.final.bam \
 bam/${base}_L_r3.t2t.NM4.sorted.final.bam bam/${base}_H_r3.t2t.NM4.sorted.final.bam


