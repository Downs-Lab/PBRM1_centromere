#!/bin/bash
#SBATCH --array=0
#SBATCH --cpus-per-task=48
#SBATCH --time=6:00:00

# Set array for number of samples of interest  e.g. #SBATCH --array=0-8 for 9 samples

# requires approx 384G (e.g. 48 threads and 8G) and 12 hours running time

cd ..

INPUT_FILES=(kmer_centromere/kmer_dbs/51/*_r1_normkcounts.txt)

base=$(basename -s _H_r1_normkcounts.txt ${INPUT_FILES[$SLURM_ARRAY_TASK_ID]})

Sample_1_count=$(cat kmer_centromere/kmer_dbs/${base}_r1_readlens.csv)
Sample_2_count=$(cat kmer_centromere/kmer_dbs/${base}_r2_readlens.csv)
Sample_3_count=$(cat kmer_centromere/kmer_dbs/${base}_r3_readlens.csv)


# Set base2 to the name of your control sample

base2="WT_IgG"

Control_1_count=$(cat kmer_centromere/kmer_dbs/${base2}_r1_readlens.csv)
Control_2_count=$(cat kmer_centromere/kmer_dbs/${base2}_r2_readlens.csv)
Control_3_count=$(cat kmer_centromere/kmer_dbs/${base2}_r3_readlens.csv)

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun Rscript k-mer_analysis_pipeline/24_enriched_kmers_3control.R $base $Sample_1_count $Sample_2_count $Sample_3_count \
 $Control_1_count $Control_2_count $Control_3_count $base2

