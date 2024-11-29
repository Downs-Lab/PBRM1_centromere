#!/bin/bash
#SBATCH --cpus-per-task=6
#SBATCH --time=24:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8000


# for loop used here instead of array - otherwise would need to set up tmp folders for each array separately - or causes problems with kmc script


# Identify k-mer sequences and counts in all reads (filtered to remove spike-in reads), and then normalise against total base count

# requires approx 48G (e.g. 6 threads and 8G) and <24 hours running time


for files in filtered_fastq/*_cen_filter_1.fq.gz;
do

base=$(basename -s _cen_filter_1.fq.gz $files)


#rid_PCR_dups:

clumpify.sh t=6 in=filtered_fastq/${base}_cen_filter_1.fq.gz in2=filtered_fastq/${base}_cen_filter_2.fq.gz \
out=kmer_centromere/deduped/${base}_cen_filter_1.fq.gz out2=kmer_centromere/deduped/${base}_cen_filter_2.fq.gz dedupe=t subs=0 \
reorder=f overwrite=t -Xmx6g 2>kmer_centromere/deduped/${base}_deduping.log

#rule trim_adapters_PE:

## download adapter .fa file from here - https://github.com/timflutre/trimmomatic/tree/master/adapters - can also try this file which has more extensive adapter options (TruSeq3-PE-2.fa)

trimmomatic -XX:ParallelGCThreads=6 -Xmx6g PE -phred33 kmer_centromere/deduped/${base}_cen_filter_1.fq.gz \
kmer_centromere/deduped/${base}_cen_filter_2.fq.gz kmer_centromere/trim/${base}_paired_1.fq.gz kmer_centromere/trim/${base}_unpaired_1.fq.gz \
kmer_centromere/trim/${base}_paired_2.fq.gz kmer_centromere/trim/${base}_unpaired_2.fq.gz ILLUMINACLIP:genomes/TruSeq3-PE.fa:2:30:12 \
SLIDINGWINDOW:10:10 MINLEN:71

#rule pair_reads:

reformat.sh in1=kmer_centromere/trim/${base}_paired_1.fq.gz in2=kmer_centromere/trim/${base}_paired_2.fq.gz \
out=kmer_centromere/interleaved/${base}_paired.fq.gz

#rule feed_kmc:

rm -rf "kmer_centromere/kmer_dbs"/tmp
mkdir -p "kmer_centromere/kmer_dbs"/tmp
ln -s $(readlink -f kmer_centromere/interleaved/${base}_paired.fq.gz) kmer_centromere/kmer_dbs/${base}.fq.gz

# rule generate_kmers ## can try various kmer lengths here:

kmc -k51 -sm -t6 -ci2 -cs100000000 kmer_centromere/kmer_dbs/${base}.fq.gz kmer_centromere/kmer_dbs/51/${base}_kcounts "kmer_centromere/kmer_dbs"/tmp

#rule compute_total_bases:

zcat kmer_centromere/interleaved/${base}_paired.fq.gz | sed -e 's/\t/ /g' | \
awk -F '[ ]' 'BEGIN{{OFS=" ";}}(NR%4==1){{x=$0; y=(NR+3)/4; getline; printf("%d\n",length($1));}}' | \
awk -F '/n' '{{sum+=$1}} END {{print sum}}' > kmer_centromere/kmer_dbs/${base}_readlens.csv

# rule dump_kmers:

kmc_tools transform kmer_centromere/kmer_dbs/51/${base}_kcounts dump -s kmer_centromere/kmer_dbs/51/${base}_kcounts.dump.sort

# rule norm_kcounts (in my case I normalised centromere kmers against kmers from fq.gz files with spike in reads removed):

mylen="$(cat kmer_all/kmer_dbs/${base}_readlens.csv)"
awk -v x=$mylen '{{print$0, $2/x}}' kmer_centromere/kmer_dbs/51/${base}_kcounts.dump.sort > kmer_centromere/kmer_dbs/51/${base}_normkcounts.txt

done
