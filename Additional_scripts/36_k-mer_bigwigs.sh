#!/bin/bash
#SBATCH --array=0-8
#SBATCH --cpus-per-task=6
#SBATCH --time=8:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8000

# array currently set for data from publication
# e.g. #SBATCH --array=0-8 for 9 samples (WT_PBRM1, WT_SMARCA4, PBRM1KO_SMARCA4, WT_CENPB, SMARCA4KO_SMARCA4, WT_H3K9me2, WT_H3K9me3,
#PBRM1KO_H3K9me2, PBRM1KO_H3K9me3)

cd ..

mkdir -p kmer_centromere/enriched_kmers/overlap_2/bigwig

INPUT_FILES=(kmer_centromere/enriched_kmers/overlap_2/bam/*.bam)

base=$(basename -s .bam ${INPUT_FILES[$SLURM_ARRAY_TASK_ID]})

mkdir -p unique_reads/bigwig

# index unique bam files

samtools index kmer_centromere/enriched_kmers/overlap_2/bam/${base}.bam

# create bigwig files using scale factor calculated above

bamCoverage -p 6 --smoothLength 1 --binSize 1 --scaleFactor 1 --outFileFormat bigwig \
-b kmer_centromere/enriched_kmers/overlap_2/bam/${base}.bam -o kmer_centromere/enriched_kmers/overlap_2/bigwig/${base}.bs10.bw

