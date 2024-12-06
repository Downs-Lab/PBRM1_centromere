#!/bin/bash
#SBATCH --array=0-8
#SBATCH --cpus-per-task=6
#SBATCH --time=4:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8000

# array currently set for data from publication
# e.g. #SBATCH --array=0-8 for 9 samples (WT_PBRM1, WT_SMARCA4, PBRM1KO_SMARCA4, WT_CENPB, SMARCA4KO_SMARCA4, WT_H3K9me2, WT_H3K9me3,
#PBRM1KO_H3K9me2, PBRM1KO_H3K9me3)


#Map enriched kmers back to genome, allowing up to 5000 multimapping locations - with no mismatches

# requires approx 48G (e.g. 6 threads and 8G) and <4 hour running time

cd ..

INDEX=genomes/human_genome_t2tv1.1_chr_indexed

INPUT_FILES=(kmer_centromere/enriched_kmers/overlap_2/*.fa)

base=$(basename -s .fa ${INPUT_FILES[$SLURM_ARRAY_TASK_ID]})

echo ${INPUT_FILES[$SLURM_ARRAY_TASK_ID]} >> kmer_centromere/enriched_kmers/overlap_2/bam/bowtie_log.txt
   bowtie2 -f -x $INDEX -p 6 -U kmer_centromere/enriched_kmers/overlap_2/${base}.fa -k 5000 --score-min C,0,0 \
   2>>kmer_centromere/enriched_kmers/overlap_2/bam/bowtie_log.txt | samtools view -S -b -o kmer_centromere/enriched_kmers/overlap_2/bam/$base.bam
   echo "sorting $base.bam" >> kmer_centromere/enriched_kmers/overlap_2/bam/bowtie_log.txt
   samtools sort kmer_centromere/enriched_kmers/overlap_2/bam/$base.bam -T $base -o kmer_centromere/enriched_kmers/overlap_2/bam/$base.sort.bam

mv kmer_centromere/enriched_kmers/overlap_2/bam/$base.sort.bam kmer_centromere/enriched_kmers/overlap_2/bam/$base.bam

