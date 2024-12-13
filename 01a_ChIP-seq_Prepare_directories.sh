#!/bin/bash

# Set up directories and prepare files for ChIP-seq datasets analysed in the publication

# set custom path to initial folder where script is

cd /path/to/file/CEN_analysis

# Prepare folder for genome files and subfolder for combined genome (prepared later)

mkdir -p genomes/combined_t2tsc

# Download following files to genomes folder (not combined_t2tsc subfolder)
# CHM13-T2T and yeast spike-in genomes, cenSatAnnotation file and Trimmomatic adapter file
# S. cerevisiae: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000146045.2/
# CHM13-T2T v1.1: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_009914755.3/
# cenSatAnnotation: http://t2t.gi.ucsc.edu/chm13/hub/t2t-chm13-v1.0/cenSatAnnotation.bigBed
# ILLUMINACLIP for Trimmomatic adapter file - https://github.com/timflutre/trimmomatic/tree/master/adapters/TruSeq3-PE.fa

mkdir -p fastq/ChIP-seq/GSE175752/cell_cycle

mkdir -p fastq/ChIP-seq/GSE152681

mkdir -p fastq/ChIP-seq/GSE71510

# Download the following fq.gz files from GEO accession GSE175752 and add to the equivalent folder 
# (all RPE.Ctrl.* ChIP-seq samples e.g. RPE.Ctrl.H3K4me1.G1.rep1 (H3K4me1/2/3, H3K9ac/me1/2/3, H3K27ac/me3, H3K36me1/2/3, input - G1/ES/LS/G2 - rep1/rep2)
# naming system e.g. H3K4me1_G1_r1.fq.gz

# combine fq.gz files across the cell cyle to one

for files in fastq/ChIP-seq/GSE175752/cell_cycle/*.fq.gz;
do
base=$(basename -s _G1_r1.fq.gz $files)

cat fastq/ChIP-seq/GSE175752/cell_cycle/${base}*_r1.fq.gz > fastq/ChIP-seq/GSE175752/${base}_r1.fq.gz
cat fastq/ChIP-seq/GSE175752/cell_cycle/${base}*_r2.fq.gz > fastq/ChIP-seq/GSE175752/${base}_r2.fq.gz

done

# Download the fq.gz files with the following SRA accession codes from GEO accession GSE152681 and add to the equivalent folder 
# SRR1203667, SRR23588908, SRR12036684, SRR23588912, SRR12036698, SRR12036715, SRR23588903, SRR12036716, SRR12036717,
# SRR12036718, SRR12036719, SRR12036720, SRR12036721, SRR12036722, SRR12036723, SRR12036724, SRR12036725, SRR12036726, SRR12036727, SRR12036728
# naming system e.g. SRR12036678_PBRM1_786O_PBRM1WT.fq.gz, SRR12036715_PBRM1_A498.fq.gz

# combine fastq files that were run separately but had the same sample information

cat fastq/ChIP-seq/GSE152681/SRR12036678_PBRM1_786O_PBRM1WT.fq.gz fastq/ChIP-seq/GSE152681/SRR23588908_PBRM1_786O_PBRM1WT.fq.gz > fastq/ChIP-seq/GSE152681/PBRM1_786O_PBRM1WT.fq.gz

cat fastq/ChIP-seq/GSE152681/SRR12036684_SMARCA4_786O_PBRM1WT.fq.gz fastq/ChIP-seq/GSE152681/SRR23588912_SMARCA4_786O_PBRM1WT.fq.gz > fastq/ChIP-seq/GSE152681/SMARCA4_786O_PBRM1WT.fq.gz

mv fastq/ChIP-seq/GSE152681/SRR12036698_Input_SE_786O_PBRM1WT.fq.gz fastq/ChIP-seq/GSE152681/Input_SE_786O_PBRM1WT.fq.gz

cat fastq/ChIP-seq/GSE152681/SRR12036715_PBRM1_A498.fq.gz fastq/ChIP-seq/GSE152681/SRR23588903_PBRM1_A498.fq.gz > fastq/ChIP-seq/GSE152681/PBRM1_A498.fq.gz

mv fastq/ChIP-seq/GSE152681/SRR12036716_Input_A498.fq.gz fastq/ChIP-seq/GSE152681/Input_A498.fq.gz

cat fastq/ChIP-seq/GSE152681/SRR12036717_PBRM1_HK2.fq.gz fastq/ChIP-seq/GSE152681/SRR12036718_PBRM1_HK2.fq.gz fastq/ChIP-seq/GSE152681/SRR12036719_PBRM1_HK2.fq.gz \
fastq/ChIP-seq/GSE152681/SRR12036720_PBRM1_HK2.fq.gz fastq/ChIP-seq/GSE152681/SRR12036721_PBRM1_HK2.fq.gz fastq/ChIP-seq/GSE152681/SRR12036722_PBRM1_HK2.fq.gz \
> fastq/ChIP-seq/GSE152681/PBRM1_HK2.fq.gz

cat fastq/ChIP-seq/GSE152681/SRR12036723_Input_HK2.fq.gz fastq/ChIP-seq/GSE152681/SRR12036724_Input_HK2.fq.gz fastq/ChIP-seq/GSE152681/SRR12036725_Input_HK2.fq.gz \
fastq/ChIP-seq/GSE152681/SRR12036726_Input_HK2.fq.gz fastq/ChIP-seq/GSE152681/SRR12036727_Input_HK2.fq.gz fastq/ChIP-seq/GSE152681/SRR12036728_Input_HK2.fq.gz \
> fastq/ChIP-seq/GSE152681/Input_HK2.fq.gz

rm fastq/ChIP-seq/GSE152681/SRR*

# Download the following fq.gz files from GEO accession GSE71510 and add to the equivalent folder 
# ChIP-seq_HCT116_Parental_SMARCA4III, ChIP-seq_HCT116_Parental_inputIII
# naming system e.g. HCT116_Parental_SMARCA4.fq.gz

# prepare other folders necessary for downstream analysis

mkdir -p fastqc/ChIP-seq/GSE175752/trim
mkdir -p fastqc/ChIP-seq/GSE152681/trim
mkdir -p fastqc/ChIP-seq/GSE71510/trim

mkdir -p trimmed_reads/ChIP-seq/GSE175752
mkdir -p trimmed_reads/ChIP-seq/GSE152681
mkdir -p trimmed_reads/ChIP-seq/GSE71510

mkdir -p sam/ChIP-seq/GSE175752
mkdir -p sam/ChIP-seq/GSE152681
mkdir -p sam/ChIP-seq/GSE71510

mkdir -p bam/ChIP-seq/GSE175752
mkdir -p bam/ChIP-seq/GSE152681
mkdir -p bam/ChIP-seq/GSE71510

mkdir -p unique_reads/ChIP-seq/GSE175752/MACS2_peaks/overlap_2/CEN_peaks
mkdir -p unique_reads/ChIP-seq/GSE152681/MACS2_peaks/overlap_2/CEN_peaks
mkdir -p unique_reads/ChIP-seq/GSE71510/MACS2_peaks/overlap_2/CEN_peaks

mkdir -p bam/ChIP-seq/sc_reads

mkdir -p filtered_fastq/ChIP-seq

mkdir -p kmer_all/ChIP-seq/deduped
mkdir -p kmer_all/ChIP-seq/trim
mkdir -p kmer_all/ChIP-seq/kmer_dbs/51

mkdir -p kmer_centromere/ChIP-seq/deduped
mkdir -p kmer_centromere/ChIP-seq/trim
mkdir -p kmer_centromere/ChIP-seq/kmer_dbs/51/control
mkdir -p kmer_centromere/ChIP-seq/enriched_kmers/overlap_2/bam
mkdir -p kmer_centromere/ChIP-seq/enriched_kmers/overlap_2/intersect_CEN
mkdir -p kmer_centromere/ChIP-seq/enriched_kmers/overlap_2/pivot_CEN