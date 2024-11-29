#!/bin/bash
#SBATCH --array=0-81
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=4000

# array currently set for data from publication
# e.g. #SBATCH --array=0-81 for 82 samples in array - 3x PBRM1, 3x SMARCA4, 3x CENPB, 3x H3K9me2, 3 x H3K9me3 and 8x IgG in parental RPE1,
# and 3x SMARCA4, 3x H3K9me2, 3x H3K9me3 and 6x IgG in PBRM1KO RPE1, 3x SMARCA4 in SMARCA4KO RPE1 - low and high salt for each

## Filter out reads with more than 3 mismatches and their corresponding mates ##

## Requires approx 4Gb of memory and 40 minutes running time

INPUT_FILES=(sam/*_t2tsc.sam)
base=$(basename -s _t2tsc.sam ${INPUT_FILES[$SLURM_ARRAY_TASK_ID]})

sambamba view --sam-input -h -t 1 -F "[NM] < 4" -f sam -o sam/${base}_t2tsc.NM4.sam sam/${base}_t2tsc.sam

samtools sort \
--output-fmt BAM \
-o bam/${base}.t2tsc.NM4.sorted.bam \
sam/${base}_t2tsc.NM4.sam

## extract list of unpaired reads, filter from bam file ##

awk '{ print $1 }' sam/${base}_t2tsc.NM4.sam > sam/${base}_t2tsc.NM4.col1.txt

sort sam/${base}_t2tsc.NM4.col1.txt | uniq -u > sam/${base}_t2tsc.NM4.col1.uniq.txt

tail -n +2 sam/${base}_t2tsc.NM4.col1.uniq.txt > sam/${base}_t2tsc.NM4.unpaired.txt

#set path

PICARD='/path/to/file/picard.jar'

java -Xmx4g -Djava.io.tmpdir=/tmp -XX:-UseCompressedClassPointers \
-jar $PICARD FilterSamReads \
I=bam/${base}.t2tsc.NM4.sorted.bam \
O=bam/${base}.t2tsc.NM4.sorted.final.bam \
READ_LIST_FILE=sam/${base}_t2tsc.NM4.unpaired.txt \
FILTER=excludeReadList
