#!/bin/bash
#SBATCH --array=0-7
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8000

# array currently set for data from publication
# e.g. #SBATCH --array=0-7 for 8 sets of peaks analysed in last step if analysing all CUT&RUN data from publication - PBRM1, SMARCA4, CENPB, H3K9me2,
#H3K9me3 in parental RPE1, and SMARCA4, H3K9me2, H3K9me3 in PBRM1KO RPE1

# requires <8G and <1 hour running time

cd ..

INPUT_FILES=(unique_reads/MACS2_peaks/overlap_2/*_overlap_peaks_merge.bed)

base=$(basename -s _overlap_peaks_merge.bed ${INPUT_FILES[$SLURM_ARRAY_TASK_ID]})

# find peaks that are at least in 2 reps in the centromere/peri-centromere

bedtools intersect -wa -u -a unique_reads/MACS2_peaks/overlap_2/${base}_overlap_peaks_merge.bed \
 -b genomes/cenSatAnnotation.bed > unique_reads/MACS2_peaks/overlap_2/CEN_peaks/${base}_overlap_cen_peaks_merge.bed

# count number of peaks of uniquely mapping reads in at least 2 reps in centromere/peri-centromere

echo ${base} >> unique_reads/MACS2_peaks/overlap_2/CEN_peaks/cen_peaks_log.txt

wc -l unique_reads/MACS2_peaks/overlap_2/CEN_peaks/${base}_overlap_cen_peaks_merge.bed \
 >> unique_reads/MACS2_peaks/overlap_2/CEN_peaks/cen_peaks_log.txt
