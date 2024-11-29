#!/bin/bash
#SBATCH --array=0-40
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8000

# array currently set for data from publication
# e.g. #SBATCH --array=0-40 for 41 different samples in array - 3x PBRM1, 3x SMARCA4, 3x CENPB, 3x H3K9me2, 3x H3K9me3 and 8x IgG in parental RPE1,
# and 3x SMARCA4, 3x H3K9me2, 3x H3K9me3 and 6x IgG in PBRM1KO RPE1, and 3x SMARCA4 in SMARCA4KO RPE1

# requires <8G and <1 hour running time

INPUT_FILES=(bam/merged_bam/*.t2t.NM4.sorted.final.bam)

base=$(basename -s .t2t.NM4.sorted.final.bam ${INPUT_FILES[$SLURM_ARRAY_TASK_ID]})


## Index the merged bam file for downstream analysis ##


samtools index bam/merged_bam/${base}.t2t.NM4.sorted.final.bam


#bam stats to count number of reads in a bam file - useful for calculating scale factors for bigwig file outputs

samtools stats bam/merged_bam/${base}.t2t.NM4.sorted.final.bam > bam/merged_bam/${base}.t2t.NM4.sorted.final.bam.stats.txt

grep ^SN bam/merged_bam/${base}.t2t.NM4.sorted.final.bam.stats.txt \
 | cut -f 2- > bam/merged_bam/${base}.t2t.NM4.sorted.final.bam.SN.stats.txt
