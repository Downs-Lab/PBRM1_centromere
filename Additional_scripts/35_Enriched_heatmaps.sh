#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8000

## Example script from data in publication - enriched heatmaps showing enrichment signal across cell lines for SMARCA4 and PBRM1 peaks


# Requires a separate conda environment with conda install of the following packages: conda-forge::r-base=4.3.1
# bioconductor-complexheatmap=2.18.0 bioconductor-enrichedheatmap=1.32.0 r-circlize=0.4.16
# bioconductor-rtracklayer=1.62.0 bioconductor-genomicranges=1.54.1 r-data.table=1.15.4

# compatible so can add to previous environment with below packages
# bioconda::bioconductor-genomeinfodb=1.38.1 bioconda::bioconductor-annotationdbi=1.64.1
# bioconda::bioconductor-chipseeker=1.38.0 bioconda::bioconductor-ensembldb=2.26.0 bioconda::bioconductor-org.hs.eg.db=3.18.0


cd ..

mkdir -p unique_reads/Enriched_heatmaps

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun Rscript Additional_scripts/35_Enriched_heatmaps.R

