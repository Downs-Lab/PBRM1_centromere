#!/bin/bash
#SBATCH --cpus-per-task=24
#SBATCH --time=28:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8000

cd ..

# Download motif_databases folder and software source code from: https://meme-suite.org/meme/doc/download.html
# Installation guide: https://meme-suite.org/meme/doc/install.html?man_type=web

# tested with MEME/5.4.1-gompi-2021b

# Also requires bedtools package (tested with version 2.29.2)

# Analysis of CENPB enriched motifs from fig3H

mkdir -p MEME_SEA_motifs

# unique peaks


bedtools getfasta -fi genomes/GCA_009914755.3_genomic_chr.fna -bed unique_reads/MACS2_peaks/overlap_2/CEN_peaks/WT_CENPB_vs_IgG_overlap_cen_peaks_merge.bed \
 -fo unique_reads/MACS2_peaks/overlap_2/CEN_peaks/WT_CENPB_vs_IgG_overlap_cen_peaks_merge.fa


sea --p unique_reads/MACS2_peaks/overlap_2/CEN_peaks/WT_CENPB_vs_IgG_overlap_cen_peaks_merge.fa \
 --m motif_databases/HOCOMOCO/H12CORE_meme_format.meme --o MEME_SEA_motifs/WT_CENPB_peaks_MEME_HUMAN_H12core_full_shuffled

# k-mers

sea --p kmer_centromere/enriched_kmers/overlap_2/WT_CENPB_overlap_2.fa \
 --m motif_databases/HOCOMOCO/H12CORE_meme_format.meme --o MEME_SEA_motifs/WT_CENPB_k-mers_MEME_HUMAN_H12core_full_shuffled

# to get CENPB motif image

meme2images motif_databases/HOCOMOCO/H12CORE_meme_format.meme meme_motif_images -motif MEME_SEA_motifs/CENPB.H12CORE.0.S.B
