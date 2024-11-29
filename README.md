
# PBRM1 directs PBAF to pericentromeres and protects centromere integrity

### Karen A. Lane, Alison Harrod, Lillian Wu, Theodoros I. Roumeliotis, Hugang Feng, Shane Foo, Katheryn A. G. Begg, Federica Schiavoni, Frank T. Zenke, Alan A. Melcher, Jyoti S. Choudhary, Jessica A. Downs
### bioRxiv 2024.06.29.601326; doi: https://doi.org/10.1101/2024.06.29.601326

Custom pipeline for analysing mapping data (e.g. CUT&RUN (.salt) or ChIP-seq) at centromeric and pericentromeric sequences - extracting uniquely mapping reads and enriched k-mer analysis 
(based on pipeline used in Smith et al. 2021, Genome research. doi:10.1101/gr.267781.120)

### 1. System requirements 

### 1.1 Software dependencies and operating systems 

### Software dependencies: 

The analysis pipeline was performed using the CentOS8 linux distribution platform (based on Redhat Enterprise Linux Operating System)
The analysis pipeline can be run on a computer with the linux operating system


### 1.2 Required non-standard hardware 

The analysis pipeline was performed on a high performance computer with 1248x Intel Xeon Platinum 8260(Cascade Lake) @ 2.40GHz cores with minimum 8GB/core
The analysis pipeline requires a minimum of 8 cores with minimum 16Gb/core 


### 2. Installation guide

Prepare necessary conda/mamba environments for running downstream script

#mamba create --name unique_reads
with following package versions: trim-galore=0.6.6 bowtie2=2.4.2 sambamba=0.7.1 picard-tools=2.23.8 SAMtools=1.11 BAMtools=2.5.1 macs2=2.2.7.1 bedtools=2.29.2 ucsc-bigbedtobed=377 FastQC=0.11.9 

#mamba create --name identify_kmers
with following package versions: Bbmap=38.87 bedtools=2.29.2 SAMtools=1.11 bowtie2=2.2.5 kmc=3.2.1 pear=0.9.6 python=3.6.13 trimmomatic=0.39 seqkit=2.5.1 datamash=1.1.0

#mamba create --name enriched_kmers
with following package versions: r-base= 4.2.3  r-dplyr=1.1.1 r-tibble= 3.2.0 python=3.11.0 r-r.utils=2.12.3 r-ggplot2=3.5.1

#mamba create --name kmer_analysis
with following package versions: bowtie2=2.4.2 SAMtools=1.11 bedtools=2.29.2 datamash=1.1.0 deeptools=3.5.1


### 3. Instructions for use 

For setup instructions refer to the file 01a_Prepare_directories.sh first, then download the files below and run 01b_Prepare_files.sh

Download following files from links into the 'genomes' folder:
CHM13-T2T and yeast spike-in genomes, cenSatAnnotation file and Trimmomatic adapter file
S. cerevisiae: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000146045.2/
CHM13-T2T v1.1: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_009914755.3/
cenSatAnnotation: http://t2t.gi.ucsc.edu/chm13/hub/t2t-chm13-v1.0/cenSatAnnotation.bigBed
ILLUMINACLIP for Trimmomatic adapter file - https://github.com/timflutre/trimmomatic/tree/master/adapters/TruSeq3-PE.fa

Download or prepare data files:
If using data from the publication, then download fq.gz files from GEO accession GSE235294 and add to the 'fastq' folder, with same naming system as metadata
e.g. WT_SMARCA4_L_r1_1.fq.gz
or add your own fq.gz files ready for analysis to this folder

Use of conda environments created:
unique_reads is used for 02_trim_fastqc.sh to 16_cen_peaks.sh
identify_kmers is used for 18_Filter_out_spike_in.sh to 23_mv_control_kmers.sh
enriched_kmers is used for 24_enriched_kmers_9IgG_LH.sh, 25_overlap_2reps_enriched_kmers.sh, 29_kmer_specific_groups.sh
kmer_analysis is used for 26_enrich_kmer_fa.sh to 28_cen_pivot.sh
 
For unique read mapping:
follow 02_trim_fastqc.sh to 17_specific_groups.sh in numerical order 
(within 'Unique_read_mapping' folder, including 7b_mark_remove_dups.sh)

For kmer analysis: 
Output files generated from steps up to 12_index_stats.sh are required
Follow 18_Filter_out_spike_in.sh to 29_kmer_specific_groups.sh (within 'k-mer_analysis_pipeline' folder)

## Citation

Please cite the following if you find this analysis pipeline useful:

PBRM1 directs PBAF to pericentromeres and protects centromere integrity
Karen A. Lane, Alison Harrod, Lillian Wu, Theodoros I. Roumeliotis, Hugang Feng, Shane Foo, Katheryn A. G. Begg, Federica Schiavoni, Frank T. Zenke, Alan A. Melcher, Jyoti S. Choudhary, Jessica A. Downs
bioRxiv 2024.06.29.601326; doi: https://doi.org/10.1101/2024.06.29.601326