#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8000

# requires the following conda package intervene=0.6.5

cd ..

# set base and base2 to sample names you want to compare

base="Sample_1_vs_IgG"

base2="Sample_2_vs_IgG"

mkdir -p unique_reads/MACS2_peaks/overlap_2/intervene

intervene venn -i unique_reads/MACS2_peaks/overlap_2/${base}_overlap_peaks_merge.bed unique_reads/MACS2_peaks/overlap_2/${base2}_overlap_peaks_merge.bed \
 --names=${base},${base2} --output unique_reads/MACS2_peaks/overlap_2/intervene/${base}_${base2}_intervene

mkdir -p unique_reads/MACS2_peaks/overlap_2/CEN_peaks/intervene

intervene venn -i unique_reads/MACS2_peaks/overlap_2/CEN_peaks/${base}_overlap_cen_peaks_merge.bed \
unique_reads/MACS2_peaks/overlap_2/CEN_peaks/${base2}_overlap_cen_peaks_merge.bed  --names=${base},${base2} \
--output unique_reads/intervene/${base}_${base2}_intervene
