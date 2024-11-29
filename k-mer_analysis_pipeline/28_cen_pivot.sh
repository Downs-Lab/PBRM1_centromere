#!/bin/bash
#SBATCH --array=0-8
#SBATCH --cpus-per-task=1
#SBATCH --time=36:00:00
#SBATCH --partition=smp
#SBATCH --mem-per-cpu=16000

# array currently set for data from publication
# e.g. #SBATCH --array=0-8 for 9 samples (WT_PBRM1, WT_SMARCA4, PBRM1KO_SMARCA4, WT_CENPB, SMARCA4KO_SMARCA4, WT_H3K9me2, WT_H3K9me3,
#PBRM1KO_H3K9me2, PBRM1KO_H3K9me3)

# requires approx 16G (unable to use multithreaded job with datamash) and <36 hour running time


INPUT_FILES=(kmer_centromere/enriched_kmers/overlap_2/bam/*.bam)

base=$(basename -s .bam ${INPUT_FILES[$SLURM_ARRAY_TASK_ID]})

bedtools intersect -abam kmer_centromere/enriched_kmers/overlap_2/bam/${base}.bam -b genomes/cenSatAnnotation.bed -bed -wb \
 > kmer_centromere/enriched_kmers/overlap_2/intersect_CEN/${base}_CEN_details.bed


# prepare pivot table to summarise number of kmers in each part of annotation
# CEN_pivot.txt removes the first "_" in the annotation description to group by i.e. hor, ct, hsat3 etc.
# CEN_details_pivot.txt - shows how many kmers are in each specific part of the annotation

## split columns based on "_" and tabs, then take column containing beginning of annotation e.g. hor,ct etc. and create pivot table

awk -F ["_","\t"] 'OFS="\t" {print $4,$20}' kmer_centromere/enriched_kmers/overlap_2/intersect_CEN/${base}_CEN_details.bed  | sort -k 1 | \
datamash -s crosstab 1,2 | sed 's~N/A~0~g'  > kmer_centromere/enriched_kmers/overlap_2/pivot_CEN/${base}_CEN_pivot.txt

sort -k 4 kmer_centromere/enriched_kmers/overlap_2/intersect_CEN/${base}_CEN_details.bed | datamash -s crosstab 4,16 | sed 's~N/A~0~g' \
 > kmer_centromere/enriched_kmers/overlap_2/pivot_CEN/${base}_CEN_details_pivot.txt

