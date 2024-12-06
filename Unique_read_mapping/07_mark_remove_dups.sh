#!/bin/bash
#SBATCH --array=0-81
#SBATCH --cpus-per-task=24
#SBATCH --time=1:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8000

# array currently set for data from publication
# e.g. #SBATCH --array=0-81 for 82 samples in array - 3x PBRM1, 3x SMARCA4, 3x CENPB, 3x H3K9me2, 3 x H3K9me3 and 8x IgG in parental RPE1, and 3x SMARCA4,
# 3x H3K9me2, 3x H3K9me3 and 6x IgG in PBRM1KO RPE1, 3x SMARCA4 in SMARCA4KO RPE1 - low and high salt for each


## To mark duplicates in a bam file ##
## Usually not needed for CUT&RUN libraries, unless duplication rate is high (>50%), then removing duplicates could improve data quality and downstream analysis ##
## Percentage duplication can be found in the metrics .txt file for making a decision ##
# remove duplicate script is just shown for specific samples which had duplicates removed from in this publication

# Requires approx 128GB memory and 30 minutes running time (24 threads)

cd ..

INPUT_FILES=(bam/*.t2t.NM4.sorted.final.bam)

base=$(basename -s .t2t.NM4.sorted.final.bam ${INPUT_FILES[$SLURM_ARRAY_TASK_ID]})

#set path

PICARD='/path/to/file/picard.jar'

java -Xmx128g -Djava.io.tmpdir=/tmp -XX:ParallelGCThreads=24 -XX:-UseCompressedClassPointers \
-XX:ConcGCThreads=24 -jar $PICARD MarkDuplicates \
I=bam/${base}.t2t.NM4.sorted.final.bam \
O=bam/${base}.t2t.NM4.sorted.final.mrkDup.bam \
REMOVE_DUPLICATES=false \
METRICS_FILE=bam/${base}.t2t.NM4.sorted.finalmrkDup.txt


