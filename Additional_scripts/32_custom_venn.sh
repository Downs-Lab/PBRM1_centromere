#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8000

## Example script in R file is from publication data for parental and PBRM1-KO RPE1 H3K9me2 peak intersection
## This was the script used to create the final venn diagrams using the numbers generated with the intervene script

cd ..

mkdir -p unique_reads/custom_venn

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun Rscript 32_custom_venn.R

