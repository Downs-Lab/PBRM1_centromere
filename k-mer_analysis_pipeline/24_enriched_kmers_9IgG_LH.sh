#!/bin/bash
#SBATCH --array=0-8
#SBATCH --cpus-per-task=48
#SBATCH --time=16:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8000

# array currently set for data from publication
# e.g. #SBATCH --array=0-8 for 9 samples (WT_PBRM1, WT_SMARCA4, PBRM1KO_SMARCA4, WT_CENPB, SMARCA4KO_SMARCA4, WT_H3K9me2, WT_H3K9me3,
#PBRM1KO_H3K9me2, PBRM1KO_H3K9me3)


# requires approx 384G (e.g. 48 threads and 8G) and 12 hours running time 

INPUT_FILES=(kmer_centromere/kmer_dbs/51/*_H_r1_normkcounts.txt)

base=$(basename -s _H_r1_normkcounts.txt ${INPUT_FILES[$SLURM_ARRAY_TASK_ID]})

Sample_L1_count=$(cat kmer_centromere/kmer_dbs/${base}_L_r1_readlens.csv)
Sample_H1_count=$(cat kmer_centromere/kmer_dbs/${base}_H_r1_readlens.csv)
Sample_L2_count=$(cat kmer_centromere/kmer_dbs/${base}_L_r2_readlens.csv)
Sample_H2_count=$(cat kmer_centromere/kmer_dbs/${base}_H_r2_readlens.csv)
Sample_L3_count=$(cat kmer_centromere/kmer_dbs/${base}_L_r3_readlens.csv)
Sample_H3_count=$(cat kmer_centromere/kmer_dbs/${base}_H_r3_readlens.csv)

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun Rscript 24_enriched_kmers_9IgG_LH.R $base $Sample_L1_count $Sample_H1_count $Sample_L2_count $Sample_H2_count $Sample_L3_count $Sample_H3_count
