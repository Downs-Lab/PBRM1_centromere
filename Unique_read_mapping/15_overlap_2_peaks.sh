#!/bin/bash
#SBATCH --array=0-7
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8000

# array currently set for data from publication
# e.g. #SBATCH --array=0-7 for 8 sets of peaks analysed in last step if analysing all CUT&RUN data from publication - PBRM1, SMARCA4, CENPB, H3K9me2,
#H3K9me3 in parental RPE1, and SMARCA4, H3K9me2, H3K9me3 in PBRM1KO RPE1

# requires 8G and <1 hour running time

INPUT_FILES=(unique_reads/MACS2_peaks/*_vs_IgG_r1_peaks.broadPeak)

base=$(basename -s _vs_IgG_r1_peaks.broadPeak ${INPUT_FILES[$SLURM_ARRAY_TASK_ID]})

# find regions of peaks that intersect in multiple reps

bedtools multiinter -i unique_reads/MACS2_peaks/${base}_vs_IgG_r1_peaks.broadPeak \
 unique_reads/MACS2_peaks/${base}_vs_IgG_r2_peaks.broadPeak unique_reads/MACS2_peaks/${base}_vs_IgG_r3_peaks.broadPeak \
 > unique_reads/MACS2_peaks/${base}_vs_IgG_multi_peaks.broadPeak

# extract regions that intersect in at least 2 reps

awk '$5 == "1,2" || $5 == "1,3" || $5 == "2,3" ||  $5 == "1,2,3"' \
 unique_reads/MACS2_peaks/${base}_vs_IgG_multi_peaks.broadPeak \
 > unique_reads/MACS2_peaks/${base}_vs_IgG_overlap_peaks.broadPeak

# extract first 3 columns of peak file

awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' unique_reads/MACS2_peaks/${base}_vs_IgG_overlap_peaks.broadPeak \
 > unique_reads/MACS2_peaks/${base}_vs_IgG_overlap_peaks.bed

# merge adjacent peaks that don't have a gap between (0bp distance)

bedtools merge -i unique_reads/MACS2_peaks/${base}_vs_IgG_overlap_peaks.bed \
 > unique_reads/MACS2_peaks/overlap_2/${base}_vs_IgG_overlap_peaks_merge.bed

# count number of peaks of uniquely mapping reads in at least 2 reps across whole genome

echo ${base} >> unique_reads/MACS2_peaks/overlap_2/genome_peaks_log.txt

wc -l unique_reads/MACS2_peaks/overlap_2/${base}_vs_IgG_overlap_peaks_merge.bed \
 >> unique_reads/MACS2_peaks/overlap_2/genome_peaks_log.txt
