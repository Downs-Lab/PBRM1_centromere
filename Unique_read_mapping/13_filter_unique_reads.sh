#!/bin/bash
#SBATCH --array=0-40
#SBATCH --cpus-per-task=6
#SBATCH --time=3:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8000

# array currently set for data from publication
# e.g. #SBATCH --array=0-40 for 41 different samples in array - 3x PBRM1, 3x SMARCA4, 3x CENPB, 3x H3K9me2, 3x H3K9me3 and 8x IgG in parental RPE1,
#and 3x SMARCA4, 3x H3K9me2, 3x H3K9me3 and 6x IgG in PBRM1KO RPE1, and 3x SMARCA4 in SMARCA4KO RPE1

# requires approx 48G (e.g. 6 threads and 8G) and 2 hours running time

cd ..

INPUT_FILES=(bam/merged_bam/*.t2t.NM4.sorted.final.bam)

base=$(basename -s .t2t.NM4.sorted.final.bam ${INPUT_FILES[$SLURM_ARRAY_TASK_ID]})

sambamba view -h -t 6 -f bam -F "[XS] == null and not unmapped and not duplicate" bam/merged_bam/${base}.t2t.NM4.sorted.final.bam \
> unique_reads/${base}.bam


