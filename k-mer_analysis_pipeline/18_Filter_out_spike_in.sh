#!/bin/bash
#SBATCH --array=0-81
#SBATCH --cpus-per-task=6
#SBATCH --time=12:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8000

# array currently set for data from publication
# e.g. #SBATCH --array=0-81 for 82 samples in array - 3x PBRM1, 3x SMARCA4, 3x CENPB, 3x H3K9me2, 3 x H3K9me3 and 8x IgG in parental RPE1, 
# and 3x SMARCA4, 3x H3K9me2, 3x H3K9me3 and 6x IgG in PBRM1KO RPE1, 3x SMARCA4 in SMARCA4KO RPE1 - low and high salt for each

# convert sc (yeast mapped) bams to fastq files - or alternative spike in genome

# requires approx 48G (e.g. 6 threads and 8G) and 12 hours running time

INPUT_FILES=(bam/*.sc.NM4.sorted.final.bam)

base=$(basename -s .bam ${INPUT_FILES[$SLURM_ARRAY_TASK_ID]})

bedtools bamtofastq -i bam/${base}.bam -fq bam/sc_reads/${base}.fq


# extract read IDs to text file for yeast-mapped reads

seqkit fx2tab bam/sc_reads/${base}.fq | awk -v OFS='\t' '{array[$1]=1} END {for (readID in array) print readID}' \
 > bam/sc_reads/${base}.readsID.txt








