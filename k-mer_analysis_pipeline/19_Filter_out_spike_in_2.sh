#!/bin/bash
#SBATCH --array=0-281
#SBATCH --cpus-per-task=6
#SBATCH --time=12:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8000

# array currently set for data from publication
# e.g. #SBATCH --array=0-81 for 82 samples in array - 3x PBRM1, 3x SMARCA4, 3x CENPB, 3x H3K9me2, 3 x H3K9me3 and 8x IgG in parental RPE1,
# and 3x SMARCA4, 3x H3K9me2, 3x H3K9me3 and 6x IgG in PBRM1KO RPE1, 3x SMARCA4 in SMARCA4KO RPE1 - low and high salt for each

# then filter fastq files to exclude yeast-mapped reads

# requires approx 48G (e.g. 6 threads and 8G) and 12 hours running time

INPUT_FILES=(bam/sc_reads/*.sc.NM4.sorted.final.readsID.txt)

base=$(basename -s .sc.NM4.sorted.final.readsID.txt ${INPUT_FILES[$SLURM_ARRAY_TASK_ID]})


filterbyname.sh in=fastq/${base}_1.fq.gz out=filtered_fastq/${base}_sc_filter_1.fq.gz names=bam/sc_reads/${base}.sc.NM4.sorted.final.readsID.txt

filterbyname.sh in=fastq/${base}_2.fq.gz out=filtered_fastq/${base}_sc_filter_2.fq.gz names=bam/sc_reads/${base}.sc.NM4.sorted.final.readsID.txt

