#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8000

## Example script from data in publication

# Download file from link below to genomes folder - GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/

# Requires a separate conda environment with conda install of the following packages: conda-forge::r-base=4.3.1
# bioconda::bioconductor-rtracklayer=1.62.0 bioconda::bioconductor-genomeinfodb=1.38.1 bioconda::bioconductor-annotationdbi=1.64.1
# bioconda::bioconductor-chipseeker=1.38.0 bioconda::bioconductor-ensembldb=2.26.0 bioconda::bioconductor-org.hs.eg.db=3.18.0


cd ..

mkdir -p unique_reads/ChIPseeker

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun Rscript Additional_scripts/34_ChIPseeker.R

