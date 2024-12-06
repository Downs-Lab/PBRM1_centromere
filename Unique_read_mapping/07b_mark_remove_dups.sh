#!/bin/bash
#SBATCH --cpus-per-task=24
#SBATCH --time=1:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8000

## To mark duplicates in a bam file ##
## Usually not needed for CUT&RUN libraries, unless duplication rate is high (>50%), then removing duplicates could improve data quality and downstream analysis ##
## Percentage duplication can be found in the metrics .txt file for making a decision ##
# example remove duplicate script is just shown for specific samples which had duplicates removed from in this publication

# Requires approx 128GB memory and 30 minutes running time (24 threads)

cd ..

PICARD='/path/to/file/picard.jar'


java -Xmx128g -Djava.io.tmpdir=/tmp -XX:ParallelGCThreads=24 -XX:-UseCompressedClassPointers \
-XX:ConcGCThreads=24 -jar $PICARD MarkDuplicates \
I=bam/WT_SMARCA4_L_r1.t2t.NM4.sorted.final.bam \
O=bam/WT_SMARCA4_L_r1.t2t.NM4.sorted.final.rmDup.bam \
REMOVE_DUPLICATES=true \
METRICS_FILE=bam/WT_SMARCA4_L_r1.t2t.NM4.sorted.final.rmDup.txt

java -Xmx128g -Djava.io.tmpdir=/tmp -XX:ParallelGCThreads=24 -XX:-UseCompressedClassPointers \
-XX:ConcGCThreads=24 -jar $PICARD MarkDuplicates \
I=bam/WT_IgG_L_r3.t2t.NM4.sorted.final.bam \
O=bam/WT_IgG_L_r3.t2t.NM4.sorted.final.rmDup.bam \
REMOVE_DUPLICATES=true \
METRICS_FILE=bam/WT_IgG_L_r3.t2t.NM4.sorted.final.rmDup.txt

