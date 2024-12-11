#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8000

## Example script from data in publication - example venn diagram script for overlap of SMARCA4 and PBRM1 k-mers


# Requires the enriched_kmers conda environment (see README)


cd ..

mkdir -p kmer_centromere/enriched_kmers/overlap_2/venns

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun Rscript Additional_scripts/37_k-mer_venns.R

