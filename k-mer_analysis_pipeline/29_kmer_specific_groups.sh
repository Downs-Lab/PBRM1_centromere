#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --time=2:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8000

# array currently set for data from publication

# This script is just for reproducing the specific groups of k-mers from the publication (requires all data from GEO)

cd ..

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun Rscript k-mer_analysis_pipeline/29_kmer_specific_groups.R
