#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8000

# currently set for data from publication

# requires <8G and <1 hour running time

# not possible with demo data - specific script used with full data in publication - accesible from GEO

bedtools intersect -a unique_reads/MACS2_peaks/overlap_2/CEN_peaks/WT_PBRM1_vs_IgG_overlap_cen_peaks_merge.bed \
 -b unique_reads/MACS2_peaks/overlap_2/CEN_peaks/WT_SMARCA4_vs_IgG_overlap_cen_peaks_merge.bed -u \
 > unique_reads/MACS2_peaks/overlap_2/CEN_peaks/WT_PBRM1_and_SMARCA4_cen_peaks_merge.bed

bedtools intersect -a unique_reads/MACS2_peaks/overlap_2/CEN_peaks/WT_PBRM1_vs_IgG_overlap_cen_peaks_merge.bed \
 -b unique_reads/MACS2_peaks/overlap_2/CEN_peaks/WT_SMARCA4_vs_IgG_overlap_cen_peaks_merge.bed -v \
 > unique_reads/MACS2_peaks/overlap_2/CEN_peaks/WT_PBRM1_not_SMARCA4_cen_peaks_merge.bed

bedtools intersect -a unique_reads/MACS2_peaks/overlap_2/CEN_peaks/WT_PBRM1_and_SMARCA4_peaks_cen_peaks_merge.bed \
 -b unique_reads/MACS2_peaks/overlap_2/CEN_peaks/PBRM1KO_SMARCA4_vs_IgG_overlap_cen_peaks_merge.bed -u \
 > unique_reads/MACS2_peaks/overlap_2/CEN_peaks/WT_PBRM1_and_SMARCA4_and_PKO_cen_peaks_merge.bed

bedtools intersect -a unique_reads/MACS2_peaks/overlap_2/CEN_peaks/WT_PBRM1_and_SMARCA4_cen_peaks_merge.bed \
 -b unique_reads/MACS2_peaks/overlap_2/CEN_peaks/PBRM1KO_SMARCA4_vs_IgG_overlap_cen_peaks_merge.bed -v \
 > unique_reads/MACS2_peaks/overlap_2/CEN_peaks/WT_PBRM1_and_SMARCA4_not_PKO_cen_peaks_merge.bed

bedtools intersect -a unique_reads/MACS2_peaks/overlap_2/CEN_peaks/WT_PBRM1_not_SMARCA4_cen_peaks_merge.bed \
 -b unique_reads/MACS2_peaks/overlap_2/CEN_peaks/PBRM1KO_SMARCA4_vs_IgG_overlap_cen_peaks_merge.bed -v \
 > unique_reads/MACS2_peaks/overlap_2/CEN_peaks/WT_PBRM1_not_SMARCA4_not_PKO_cen_peaks_merge.bed

bedtools intersect -a unique_reads/MACS2_peaks/overlap_2/CEN_peaks/WT_PBRM1_vs_IgG_overlap_cen_peaks_merge.bed \
 -b unique_reads/MACS2_peaks/overlap_2/CEN_peaks/PBRM1KO_SMARCA4_vs_IgG_overlap_cen_peaks_merge.bed -u \
 > unique_reads/MACS2_peaks/overlap_2/CEN_peaks/WT_PBRM1_and_PKO_cen_peaks_merge.bed

bedtools intersect -a unique_reads/MACS2_peaks/overlap_2/CEN_peaks/WT_PBRM1_and_PKO_cen_peaks_merge.bed \
 -b unique_reads/MACS2_peaks/overlap_2/CEN_peaks/WT_SMARCA4_vs_IgG_overlap_cen_peaks_merge.bed -v \
 > unique_reads/MACS2_peaks/overlap_2/CEN_peaks/WT_PBRM1_and_PKO_not_SMARCA4_cen_peaks_merge.bed

bedtools intersect -a unique_reads/MACS2_peaks/overlap_2/CEN_peaks/PBRM1KO_SMARCA4_vs_IgG_overlap_cen_peaks_merge.bed \
 -b unique_reads/MACS2_peaks/overlap_2/CEN_peaks/WT_SMARCA4_vs_IgG_overlap_cen_peaks_merge.bed -v \
 > unique_reads/MACS2_peaks/overlap_2/CEN_peaks/PKO_not_SMARCA4_cen_peaks_merge.bed

bedtools intersect -a unique_reads/MACS2_peaks/overlap_2/CEN_peaks/PKO_not_SMARCA4_cen_peaks_merge.bed \
 -b unique_reads/MACS2_peaks/overlap_2/CEN_peaks/WT_PBRM1_vs_IgG_overlap_cen_peaks_merge.bed -v | sort -k1,1 -k2,2n \
 > unique_reads/MACS2_peaks/overlap_2/CEN_peaks/KO_specific_cen_peaks_merge.bed

cat unique_reads/MACS2_peaks/overlap_2/CEN_peaks/WT_PBRM1_not_SMARCA4_not_PKO_cen_peaks_merge.bed \
 unique_reads/MACS2_peaks/overlap_2/CEN_peaks/WT_PBRM1_and_SMARCA4_not_PKO_cen_peaks_merge.bed | sort -k1,1 -k2,2n | \
 mergeBed -i stdin > unique_reads/MACS2_peaks/overlap_2/CEN_peaks/PBRM1_dependent_cen_peaks_merge.bed

cat unique_reads/MACS2_peaks/overlap_2/CEN_peaks/WT_PBRM1_and_SMARCA4_and_PKO_cen_peaks_merge.bed \
 unique_reads/MACS2_peaks/overlap_2/CEN_peaks/WT_PBRM1_and_PKO_not_SMARCA4_cen_peaks_merge.bed | sort -k1,1 -k2,2n | \
 mergeBed -i stdin > unique_reads/MACS2_peaks/overlap_2/CEN_peaks/PBRM1_independent_cen_peaks_merge.bed

