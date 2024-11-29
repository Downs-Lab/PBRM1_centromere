#!/bin/bash
#SBATCH --array=0-81
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=4000

# array currently set for data from publication
# e.g. #SBATCH --array=0-81 for 82 samples in array - 3x PBRM1, 3x SMARCA4, 3x CENPB, 3x H3K9me2, 3 x H3K9me3 and 8x IgG in parental RPE1, 
# and 3x SMARCA4, 3x H3K9me2, 3x H3K9me3 and 6x IgG in PBRM1KO RPE1, 3x SMARCA4 in SMARCA4KO RPE1 - low and high salt for each


## Split the filtered bam file into 2 - for reads aligning to the human and yeast genomes ##
## This splits the bam file into separate files for each chromosome ##
# Requires approx 10MB memory and 15 minutes running time

INPUT_FILES=(bam/*.t2tsc.NM4.sorted.final.bam)
base=$(basename -s .t2tsc.NM4.sorted.final.bam ${INPUT_FILES[$SLURM_ARRAY_TASK_ID]})

bamtools split -in bam/${base}.t2tsc.NM4.sorted.final.bam -reference

