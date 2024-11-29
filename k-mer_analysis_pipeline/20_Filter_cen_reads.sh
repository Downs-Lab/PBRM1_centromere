#!/bin/bash
#SBATCH --array=0-81
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8000

# array currently set for data from publication
# e.g. #SBATCH --array=0-81 for 82 samples in array - 3x PBRM1, 3x SMARCA4, 3x CENPB, 3x H3K9me2, 3 x H3K9me3 and 8x IgG in parental RPE1, 
# and 3x SMARCA4, 3x H3K9me2, 3x H3K9me3 and 6x IgG in PBRM1KO RPE1, 3x SMARCA4 in SMARCA4KO RPE1 - low and high salt for each

# filter fq.gz files for reads that originally mapped to centromere and pericentromeric regions

INPUT_FILES=(bam/*.t2t.NM4.sorted.final.bam)

base=$(basename -s .t2t.NM4.sorted.final.bam ${INPUT_FILES[$SLURM_ARRAY_TASK_ID]})


samtools sort -n bam/${base}.t2t.NM4.sorted.final.bam | bedtools bamtobed -i stdin | cut -f 1,2,3,4 > filtered_fastq/$base.bed

bedtools intersect -a filtered_fastq/$base.bed -b genomes/cenSatAnnotation.bed  > filtered_fastq/$base.cen.bed

awk {'print $4'} filtered_fastq/$base.cen.bed > filtered_fastq/$base.cen_2.txt

sed 's/..$//' < filtered_fastq/$base.cen_2.txt > filtered_fastq/$base.cen.txt

filterbyname.sh in=fastq/${base}_1.fq.gz out=filtered_fastq/${base}_cen_filter_1.fq.gz names=filtered_fastq/$base.cen.txt include=t

filterbyname.sh in=fastq/${base}_2.fq.gz out=filtered_fastq/${base}_cen_filter_2.fq.gz names=filtered_fastq/$base.cen.txt include=t


