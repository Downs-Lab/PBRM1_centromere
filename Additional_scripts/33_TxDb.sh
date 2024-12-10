#!/bin/bash
#SBATCH --cpus-per-task=12
#SBATCH --time=08:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8000

## Example script from data in publication

# Download file from link below to genomes folder - GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/

# TxDb script from https://github.com/Bioconductor/txdbmaker/issues/1

# Requires a separate conda environment with conda install of the following package: conda-forge::r-base=4.4.0

cd ..

mkdir -p unique_reads/ChIPseeker

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun Rscript Additional_scripts/33_TxDb.R

