#!/bin/bash
#SBATCH --array=0-8
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8000

# array currently set for data from publication
# e.g. #SBATCH --array=0-8 for 9 samples (WT_PBRM1, WT_SMARCA4, PBRM1KO_SMARCA4, WT_CENPB, SMARCA4KO_SMARCA4, WT_H3K9me2, WT_H3K9me3,
#PBRM1KO_H3K9me2, PBRM1KO_H3K9me3)

# requires <8G and <1 hours running time

INPUT_FILES=(kmer_centromere/enriched_kmers/overlap_2/*_overlap_2.txt)

base=$(basename -s _overlap_2.txt ${INPUT_FILES[$SLURM_ARRAY_TASK_ID]})

sed '1d' kmer_centromere/enriched_kmers/overlap_2/${base}_overlap_2.txt > kmer_centromere/enriched_kmers/overlap_2/${base}_overlap_2.fa.txt

awk -F '[ ]' 'BEGIN{{OFS="\n"}}{{n=NR; x=">"n; print x, $1}}' kmer_centromere/enriched_kmers/overlap_2/${base}_overlap_2.fa.txt \
> kmer_centromere/enriched_kmers/overlap_2/${base}_overlap_2.fa
