#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8000

cd ..


# Also requires bedtools package (tested with version 2.29.2)

# Specific relative distance calculations used to prepare supplementary figure S10C-D
# _summary.txt file outputs were used to calculate fractions (count/total) and cumulative fractions for each relative distance (reldist) and then plot in Graphpad prism

# This requires rerunning 26_enrich_kmer_fa and 27_map_kmers scripts but on the output files from 29_kmer_specific_groups script

mkdir -p kmer_centromere/enriched_kmers/overlap_2/reldist

# convert bam files to bed file

bedtools bamtobed -i kmer_centromere/enriched_kmers/overlap_2/PBRM1_specific_kmers_fc2.bam > kmer_centromere/enriched_kmers/overlap_2/PBRM1_specific_kmers_fc2.bed

bedtools bamtobed -i kmer_centromere/enriched_kmers/overlap_2/PBRM1_non-specific_kmers_fc2.bam > kmer_centromere/enriched_kmers/overlap_2/PBRM1_non-specific_kmers_fc2.bed

bedtools bamtobed -i kmer_centromere/enriched_kmers/overlap_2/KO_specific_kmers_fc2.bam > kmer_centromere/enriched_kmers/overlap_2/KO_specific_kmers_fc2.bed

# Do the same as above with equivalent for histone PTM ChIP-seq GSE175752
# (after following full k-mer analysis pipeline as before, but with exceptions found in Additional_scripts/Exceptions_for_ChIP-seq.txt)

bedtools reldist -a kmer_centromere/enriched_kmers/overlap_2/PBRM1_specific_kmers_fc2.bed \
 -b kmer_centromere/ChIP_seq/enriched_kmers/overlap_2/kmers_H3K4me3_fc2.bed \
 > kmer_centromere/enriched_kmers/overlap_2/reldist/PBRM1_specific_reldist_kmers_H3K4me3_summary.txt

bedtools reldist -a kmer_centromere/enriched_kmers/overlap_2/PBRM1_non-specific_kmers_fc2.bed \
 -b kmer_centromere/ChIP_seq/enriched_kmers/overlap_2/kmers_H3K4me3_fc2.bed \
 > kmer_centromere/enriched_kmers/overlap_2/reldist/PBRM1_non-specific_reldist_kmers_H3K4me3_summary.txt

bedtools reldist -a kmer_centromere/enriched_kmers/overlap_2/KO_specific_kmers_fc2.bed \
 -b kmer_centromere/ChIP_seq/enriched_kmers/overlap_2/kmers_H3K4me3_fc2.bed \
 > kmer_centromere/enriched_kmers/overlap_2/reldist/KO_specific_reldist_kmers_H3K4me3_summary.txt


bedtools reldist -a kmer_centromere/enriched_kmers/overlap_2/PBRM1_specific_kmers_fc2.bed \
 -b kmer_centromere/ChIP_seq/enriched_kmers/overlap_2/kmers_H3K9me3_fc2.bed \
 > kmer_centromere/enriched_kmers/overlap_2/reldist/PBRM1_specific_reldist_kmers_H3K9me3_summary.txt

bedtools reldist -a kmer_centromere/enriched_kmers/overlap_2/PBRM1_non-specific_kmers_fc2.bed \
 -b kmer_centromere/ChIP_seq/enriched_kmers/overlap_2/kmers_H3K9me3_fc2.bed \
 > kmer_centromere/enriched_kmers/overlap_2/reldist/PBRM1_non-specific_reldist_kmers_H3K9me3_summary.txt

bedtools reldist -a kmer_centromere/enriched_kmers/overlap_2/KO_specific_kmers_fc2.bed \
 -b kmer_centromere/ChIP_seq/enriched_kmers/overlap_2/kmers_H3K9me3_fc2.bed \
 > kmer_centromere/enriched_kmers/overlap_2/reldist/KO_specific_reldist_kmers_H3K9me3_summary.txt


bedtools reldist -a kmer_centromere/enriched_kmers/overlap_2/PBRM1_specific_kmers_fc2.bed \
 -b kmer_centromere/ChIP_seq/enriched_kmers/overlap_2/kmers_H3K9ac_fc2.bed \
 > kmer_centromere/enriched_kmers/overlap_2/reldist/PBRM1_specific_reldist_kmers_H3K9ac_summary.txt

bedtools reldist -a kmer_centromere/enriched_kmers/overlap_2/PBRM1_non-specific_kmers_fc2.bed \
 -b kmer_centromere/ChIP_seq/enriched_kmers/overlap_2/kmers_H3K9ac_fc2.bed \
 > kmer_centromere/enriched_kmers/overlap_2/reldist/PBRM1_non-specific_reldist_kmers_H3K9ac_summary.txt

bedtools reldist -a kmer_centromere/enriched_kmers/overlap_2/KO_specific_kmers_fc2.bed \
 -b kmer_centromere/ChIP_seq/enriched_kmers/overlap_2/kmers_H3K9ac_fc2.bed \
 > kmer_centromere/enriched_kmers/overlap_2/reldist/KO_specific_reldist_kmers_H3K9ac_summary.txt
