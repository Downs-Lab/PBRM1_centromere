#!/bin/bash

# Set up directories and prepare files

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

mkdir -p fastq

# Download fq.gz files from GEO accession GSE235294 and add to this folder, with same naming system as metadata
# e.g. WT_SMARCA4_L_r1_1.fq.gz
# or for now download demo data (randomised fq.gz files with same naming system as data on GEO)

# prepare other folders necessary for downstream analysis

mkdir -p fastqc/trim

mkdir -p trimmed_reads

mkdir -p sam

mkdir -p bam/merged_bam

mkdir -p bam/sc_reads

mkdir -p unique_reads/MACS2_peaks/overlap_2/CEN_peaks

mkdir -p filtered_fastq

mkdir -p kmer_all/deduped
mkdir -p kmer_all/trim
mkdir -p kmer_all/interleaved
mkdir -p kmer_all/kmer_dbs/51

mkdir -p kmer_centromere/deduped
mkdir -p kmer_centromere/trim
mkdir -p kmer_centromere/interleaved
mkdir -p kmer_centromere/kmer_dbs/51/control
mkdir -p kmer_centromere/enriched_kmers/overlap_2/bam
mkdir -p kmer_centromere/enriched_kmers/overlap_2/intersect_CEN
mkdir -p kmer_centromere/enriched_kmers/overlap_2/pivot_CEN


