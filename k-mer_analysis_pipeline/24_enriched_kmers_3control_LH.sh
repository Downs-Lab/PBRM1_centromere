#!/bin/bash
#SBATCH --array=0
#SBATCH --cpus-per-task=48
#SBATCH --time=6:00:00

# Set array to number of samples of interest e.g. SBATCH --array=0-8 for 9 samples

# requires approx 384G (e.g. 48 threads and 8G) and 12 hours running time 

cd ..

INPUT_FILES=(kmer_centromere/kmer_dbs/51/*_H_r1_normkcounts.txt)

base=$(basename -s _H_r1_normkcounts.txt ${INPUT_FILES[$SLURM_ARRAY_TASK_ID]})

Sample_L1_count=$(cat kmer_centromere/kmer_dbs/${base}_L_r1_readlens.csv)
Sample_H1_count=$(cat kmer_centromere/kmer_dbs/${base}_H_r1_readlens.csv)
Sample_L2_count=$(cat kmer_centromere/kmer_dbs/${base}_L_r2_readlens.csv)
Sample_H2_count=$(cat kmer_centromere/kmer_dbs/${base}_H_r2_readlens.csv)
Sample_L3_count=$(cat kmer_centromere/kmer_dbs/${base}_L_r3_readlens.csv)
Sample_H3_count=$(cat kmer_centromere/kmer_dbs/${base}_H_r3_readlens.csv)


# Set base2 to the name of your control sample

base2="WT_IgG"

Control_L1_count=$(cat kmer_centromere/kmer_dbs/${base2}_L_r1_readlens.csv)
Control_H1_count=$(cat kmer_centromere/kmer_dbs/${base2}_H_r1_readlens.csv)
Control_L2_count=$(cat kmer_centromere/kmer_dbs/${base2}_L_r2_readlens.csv)
Control_H2_count=$(cat kmer_centromere/kmer_dbs/${base2}_H_r2_readlens.csv)
Control_L3_count=$(cat kmer_centromere/kmer_dbs/${base2}_L_r3_readlens.csv)
Control_H3_count=$(cat kmer_centromere/kmer_dbs/${base2}_H_r3_readlens.csv)

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun Rscript k-mer_analysis_pipeline/24_enriched_kmers_3control_LH.R $base $Sample_L1_count $Sample_H1_count $Sample_L2_count $Sample_H2_count $Sample_L3_count $Sample_H3_count \
 $Control_L1_count $Control_H1_count $Control_L2_count $Control_H2_count $Control_L3_count $Control_H3_count $base2

