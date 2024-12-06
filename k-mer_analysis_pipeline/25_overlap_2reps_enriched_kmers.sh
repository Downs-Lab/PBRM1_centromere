#!/bin/bash
#SBATCH --array=0-8
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8000

# array currently set for data from publication
# e.g. #SBATCH --array=0-8 for 9 samples (WT_PBRM1, WT_SMARCA4, PBRM1KO_SMARCA4, WT_CENPB, SMARCA4KO_SMARCA4, WT_H3K9me2, WT_H3K9me3, PBRM1KO_H3K9me2, PBRM1KO_H3K9me3)

# requires approx 64G (e.g. 8 threads and 8G) and <12 hours running time

cd ..

INPUT_FILES=(kmer_centromere/enriched_kmers/*_1_vs_Control_fc2.txt)

base=$(basename -s _1_vs_Control_fc2.txt ${INPUT_FILES[$SLURM_ARRAY_TASK_ID]})

sed '1d' kmer_centromere/enriched_kmers/${base}_1_vs_Control_fc2.txt > kmer_centromere/enriched_kmers/${base}_1_vs_Control_fc2.fa.txt
sed '1d' kmer_centromere/enriched_kmers/${base}_2_vs_Control_fc2.txt > kmer_centromere/enriched_kmers/${base}_2_vs_Control_fc2.fa.txt
sed '1d' kmer_centromere/enriched_kmers/${base}_3_vs_Control_fc2.txt > kmer_centromere/enriched_kmers/${base}_3_vs_Control_fc2.fa.txt


export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun Rscript 25_overlap_2reps_enriched_kmers.R $base
