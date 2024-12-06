#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=8000

# currently set for data from publication

## requires approx 8 hours and 8GB

cd ..

# WT_SMARCA4 peaks

macs2 callpeak -t unique_reads/WT_SMARCA4_LrmDupH_r1.bam -c unique_reads/WT_IgG_LH_r1.bam -g 3054832041 -f BAMPE \
-q 0.01 -n WT_SMARCA4_vs_IgG_r1 --broad --bdg  --broad-cutoff 0.01 --keep-dup all --outdir unique_reads/MACS2_peaks

macs2 callpeak -t unique_reads/WT_SMARCA4_LH_r2.bam -c unique_reads/WT_IgG_LH_r2.bam -g 3054832041 -f BAMPE \
-q 0.01 -n WT_SMARCA4_vs_IgG_r2 --broad --bdg  --broad-cutoff 0.01 --keep-dup all --outdir unique_reads/MACS2_peaks

macs2 callpeak -t unique_reads/WT_SMARCA4_LH_r3.bam -c unique_reads/WT_IgG_LrmDupH_r3.bam -g 3054832041 -f BAMPE \
-q 0.01 -n WT_SMARCA4_vs_IgG_r3 --broad --bdg  --broad-cutoff 0.01 --keep-dup all --outdir unique_reads/MACS2_peaks

# remaining samples analysis script below, but not included in demo

# WT_CENPB peaks

#macs2 callpeak -t unique_reads/WT_CENPB_LH_r1.bam -c unique_reads/WT_IgG_LH_r1.bam -g 3054832041 -f BAMPE \
#-q 0.01 -n WT_CENPB_vs_IgG_r1 --broad --bdg  --broad-cutoff 0.01 --keep-dup all --outdir unique_reads/MACS2_peaks

#macs2 callpeak -t unique_reads/WT_CENPB_LH_r2.bam -c unique_reads/WT_IgG_LH_r2.bam -g 3054832041 -f BAMPE \
#-q 0.01 -n WT_CENPB_vs_IgG_r2 --broad --bdg  --broad-cutoff 0.01 --keep-dup all --outdir unique_reads/MACS2_peaks

#macs2 callpeak -t unique_reads/WT_CENPB_LH_r3.bam -c unique_reads/WT_IgG_LrmDupH_r3.bam -g 3054832041 -f BAMPE \
#-q 0.01 -n WT_CENPB_vs_IgG_r3 --broad --bdg  --broad-cutoff 0.01 --keep-dup all --outdir unique_reads/MACS2_peaks


# PBRM1KO_SMARCA4 peaks

#macs2 callpeak -t unique_reads/PBRM1KO_SMARCA4_LH_r1.bam -c unique_reads/PBRM1KO_IgG_LH_r1.bam -g 3054832041 -f BAMPE \
#-q 0.01 -n PBRM1KO_SMARCA4_vs_IgG_r1 --broad --bdg  --broad-cutoff 0.01 --keep-dup all --outdir unique_reads/MACS2_peaks

#macs2 callpeak -t unique_reads/PBRM1KO_SMARCA4_LH_r2.bam -c unique_reads/PBRM1KO_IgG_LH_r2.bam -g 3054832041 -f BAMPE \
#-q 0.01 -n PBRM1KO_SMARCA4_vs_IgG_r2 --broad --bdg  --broad-cutoff 0.01 --keep-dup all --outdir unique_reads/MACS2_peaks

#macs2 callpeak -t unique_reads/PBRM1KO_SMARCA4_LH_r3.bam -c unique_reads/PBRM1KO_IgG_LH_r3.bam -g 3054832041 -f BAMPE \
#-q 0.01 -n PBRM1KO_SMARCA4_vs_IgG_r3 --broad --bdg  --broad-cutoff 0.01 --keep-dup all --outdir unique_reads/MACS2_peaks

# WT_PBRM1 peaks

#macs2 callpeak -t unique_reads/WT_PBRM1_LH_r1.bam -c unique_reads/WT_IgG_LH_r4.bam -g 3054832041 -f BAMPE \
#-q 0.01 -n WT_PBRM1_vs_IgG_r1 --broad --bdg  --broad-cutoff 0.01 --keep-dup all --outdir unique_reads/MACS2_peaks

#macs2 callpeak -t unique_reads/WT_PBRM1_LH_r2.bam -c unique_reads/WT_IgG_LH_r5.bam -g 3054832041 -f BAMPE \
#-q 0.01 -n WT_PBRM1_vs_IgG_r2 --broad --bdg  --broad-cutoff 0.01 --keep-dup all --outdir unique_reads/MACS2_peaks

#macs2 callpeak -t unique_reads/WT_PBRM1_LH_r3.bam -c unique_reads/WT_IgG_LH_r6.bam -g 3054832041 -f BAMPE \
#-q 0.01 -n WT_PBRM1_vs_IgG_r3 --broad --bdg  --broad-cutoff 0.01 --keep-dup all --outdir unique_reads/MACS2_peaks

# WT_H3K9me2 peaks

#macs2 callpeak -t unique_reads/WT_H3K9me2_LH_r1.bam -c unique_reads/WT_IgG_LH_r7.bam -g 3054832041 -f BAMPE \
#-q 0.01 -n WT_H3K9me2_vs_IgG_r1 --broad --bdg  --broad-cutoff 0.01 --keep-dup all --outdir unique_reads/MACS2_peaks

#macs2 callpeak -t unique_reads/WT_H3K9me2_LH_r2.bam -c unique_reads/WT_IgG_LH_r7.bam -g 3054832041 -f BAMPE \
#-q 0.01 -n WT_H3K9me2_vs_IgG_r2 --broad --bdg  --broad-cutoff 0.01 --keep-dup all --outdir unique_reads/MACS2_peaks

#macs2 callpeak -t unique_reads/WT_H3K9me2_LH_r3.bam -c unique_reads/WT_IgG_LH_r8.bam -g 3054832041 -f BAMPE \
#-q 0.01 -n WT_H3K9me2_vs_IgG_r3 --broad --bdg  --broad-cutoff 0.01 --keep-dup all --outdir unique_reads/MACS2_peaks

# WT_H3K9me3 peaks

#macs2 callpeak -t unique_reads/WT_H3K9me3_LH_r1.bam -c unique_reads/WT_IgG_LH_r7.bam -g 3054832041 -f BAMPE \
#-q 0.01 -n WT_H3K9me3_vs_IgG_r1 --broad --bdg  --broad-cutoff 0.01 --keep-dup all --outdir unique_reads/MACS2_peaks

#macs2 callpeak -t unique_reads/WT_H3K9me3_LH_r2.bam -c unique_reads/WT_IgG_LH_r7.bam -g 3054832041 -f BAMPE \
#-q 0.01 -n WT_H3K9me3_vs_IgG_r2 --broad --bdg  --broad-cutoff 0.01 --keep-dup all --outdir unique_reads/MACS2_peaks

#macs2 callpeak -t unique_reads/WT_H3K9me3_LH_r3.bam -c unique_reads/WT_IgG_LH_r8.bam -g 3054832041 -f BAMPE \
#-q 0.01 -n WT_H3K9me3_vs_IgG_r3 --broad --bdg  --broad-cutoff 0.01 --keep-dup all --outdir unique_reads/MACS2_peaks

# PBRM1KO_H3K9me2 peaks

#macs2 callpeak -t unique_reads/PBRM1KO_H3K9me2_LH_r1.bam -c unique_reads/PBRM1KO_IgG_LH_r4.bam -g 3054832041 -f BAMPE \
#-q 0.01 -n PBRM1KO_H3K9me2_vs_IgG_r1 --broad --bdg  --broad-cutoff 0.01 --keep-dup all --outdir unique_reads/MACS2_peaks

#macs2 callpeak -t unique_reads/PBRM1KO_H3K9me2_LH_r2.bam -c unique_reads/PBRM1KO_IgG_LH_r5.bam -g 3054832041 -f BAMPE \
#-q 0.01 -n PBRM1KO_H3K9me2_vs_IgG_r2 --broad --bdg  --broad-cutoff 0.01 --keep-dup all --outdir unique_reads/MACS2_peaks

#macs2 callpeak -t unique_reads/PBRM1KO_H3K9me2_LH_r3.bam -c unique_reads/PBRM1KO_IgG_LH_r6.bam -g 3054832041 -f BAMPE \
#-q 0.01 -n PBRM1KO_H3K9me2_vs_IgG_r3 --broad --bdg  --broad-cutoff 0.01 --keep-dup all --outdir unique_reads/MACS2_peaks

# PBRM1KO_H3K9me3 peaks

#macs2 callpeak -t unique_reads/PBRM1KO_H3K9me3_LH_r1.bam -c unique_reads/PBRM1KO_IgG_LH_r4.bam -g 3054832041 -f BAMPE \
#-q 0.01 -n PBRM1KO_H3K9me3_vs_IgG_r1 --broad --bdg  --broad-cutoff 0.01 --keep-dup all --outdir unique_reads/MACS2_peaks

#macs2 callpeak -t unique_reads/PBRM1KO_H3K9me3_LH_r2.bam -c unique_reads/PBRM1KO_IgG_LH_r5.bam -g 3054832041 -f BAMPE \
#-q 0.01 -n PBRM1KO_H3K9me3_vs_IgG_r2 --broad --bdg  --broad-cutoff 0.01 --keep-dup all --outdir unique_reads/MACS2_peaks

#macs2 callpeak -t unique_reads/PBRM1KO_H3K9me3_LH_r3.bam -c unique_reads/PBRM1KO_IgG_LH_r6.bam -g 3054832041 -f BAMPE \
#-q 0.01 -n PBRM1KO_H3K9me3_vs_IgG_r3 --broad --bdg  --broad-cutoff 0.01 --keep-dup all --outdir unique_reads/MACS2_peaks
