#!/bin/bash

# prepare genome files ready for downstream analysis

sed \
's/CP068277.2/chr1/g;s/CP068276.2/chr2/g;s/CP068275.2/chr3/g;s/CP068274.2/chr4/g;s/CP068273.2/chr5/g;s/CP068272.2\
/chr6/g;s/CP068271.2/chr7/g;s/CP068270.2/chr8/g;s/CP068269.2/chr9/g;s/CP068268.2/chr10/g;s/CP068267.2/chr11/g;s/CP068266.2\
/chr12/g;s/CP068265.2/chr13/g;s/CP068264.2/chr14/g;s/CP068263.2/chr15/g;s/CP068262.2/chr16/g;s/CP068261.2/chr17/g;s/CP068260.2\
/chr18/g;s/CP068259.2/chr19/g;s/CP068258.2/chr20/g;s/CP068257.2/chr21/g;s/CP068256.2/chr22/g;s/CP068255.2/chrX/g;s/CP068254.1\
/chrM/g' genomes/GCA_009914755.3_genomic.fna > genomes/GCA_009914755.3_genomic_chr.fna

# combined the two reference genomes into one

cat genomes/GCA_009914755.3_genomic_chr.fna genomes/GCF_000146045.2_R64_genomic.fna > genomes/combined_t2tsc/CHM13-T2T_yeast.fna

# UCSC package: bigBedToBed v377 used to convert bigBed file to bed file for downstream analysis

bigBedToBed genomes/cenSatAnnotation.bigBed genomes/cenSatAnnotation.bed

# Bowtie2 build

mamba activate unique_reads

bowtie2-build genomes/combined_t2tsc/CHM13-T2T_yeast combined
bowtie2-build genomes/GCA_009914755.3_genomic_chr human_genome_t2tv1.1_chr_indexed

## Note: set path for picard.jar file (when installing picard-tools=2.23.8) in the following script files: 04_filter_bam.sh ; 07_mark_remove_dups.sh
