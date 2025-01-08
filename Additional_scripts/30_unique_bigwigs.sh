#!/bin/bash
#SBATCH --array=0-40
#SBATCH --cpus-per-task=12
#SBATCH --time=8:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8000

# array currently set for data from publication
# e.g. #SBATCH --array=0-40 for 41 different samples in array - 3x PBRM1, 3x SMARCA4, 3x CENPB, 3x H3K9me2, 3x H3K9me3 and 8x IgG in parental RPE1,
# and 3x SMARCA4, 3x H3K9me2, 3x H3K9me3 and 6x IgG in PBRM1KO RPE1, and 3x SMARCA4 in SMARCA4KO RPE1

# requires deeptools=3.5.1 - can use conda environment kmer_analysis (see README)

cd ..

INPUT_FILES=(unique_reads/*.bam)

base=$(basename -s .bam ${INPUT_FILES[$SLURM_ARRAY_TASK_ID]})

# choose scaling_reads value so that your final scaling factors are <=1

scaling_reads=10000000

# extract number of reads from bam stats file (before filtering for unique reads)

fragments=$(grep '1st fragments:' bam/merged_bam/${base}.t2t.NM4.sorted.final.bam.SN.stats.txt | cut -c 16-)

# set scale factor

scalefactor=$(bc <<<"scale=8; $scaling_reads / $fragments")

echo "number of fragments is" $fragments

echo "scale factor is" $scalefactor

mkdir -p unique_reads/bigwig

# index unique bam files

samtools index unique_reads/${base}.bam

# create bigwig files using scale factor calculated above

bamCoverage -p 12 --maxFragmentLength 1500 --extendReads --scaleFactor $scalefactor --binSize 10 \
-b unique_reads/${base}.bam -o unique_reads/bigwig/${base}.bs10.bw

