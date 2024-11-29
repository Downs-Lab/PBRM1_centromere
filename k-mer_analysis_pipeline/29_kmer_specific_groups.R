
setwd("kmer_centromere/enriched_kmers/overlap_2")

# load files

WT_SMARCA4_2 <- read.delim(file ="WT_SMARCA4_overlap_2.txt", header = F, sep = "")
WT_SMARCA4_2 <- WT_SMARCA4_2[,1]

PBRM1KO_SMARCA4_2 <- read.delim(file ="PBRM1KO_SMARCA4_overlap_2.txt", header = F, sep = "")
PBRM1KO_SMARCA4_2 <- PBRM1KO_SMARCA4_2[,1]

WT_PBRM1_2 <- read.delim(file ="WT_PBRM1_overlap_2.txt", header = F, sep = "")
WT_PBRM1_2 <- WT_PBRM1_2[,1]


# SMARCA4, PBRM1 and PBRM1-KO SMARCA4 k-mer overlap used to calculate specific groups

overlap_SMARCA4_and_PBRM1_and_PBRM1KO_SMARCA4_2 <- list(WT_SMARCA4_2 = WT_SMARCA4_2, WT_PBRM1_2 = WT_PBRM1_2, PBRM1KO_SMARCA4_2 = PBRM1KO_SMARCA4_2)

overlap_SMARCA4_and_PBRM1_and_PBRM1KO_SMARCA4_2_overlap <- overlap(Venn(overlap_SMARCA4_and_PBRM1_and_PBRM1KO_SMARCA4_2))

Just_WT_SMARCA4 <- setdiff(WT_SMARCA4_2, c(WT_PBRM1_2, PBRM1KO_SMARCA4_2))
Just_PBRM1KO_SMARCA4 <- setdiff(PBRM1KO_SMARCA4_2, c(WT_PBRM1_2, WT_SMARCA4_2))
Just_WT_PBRM1 <- setdiff(WT_PBRM1_2, c(PBRM1KO_SMARCA4_2, WT_SMARCA4_2))


overlap_SMARCA4_and_PBRM1_2 <- list(WT_SMARCA4_2 = WT_SMARCA4_2, WT_PBRM1_2 = WT_PBRM1_2)
overlap_SMARCA4_and_PBRM1_2_overlap <- overlap(Venn(overlap_SMARCA4_and_PBRM1_2))
overlap_PBRM1KO_SMARCA4_and_PBRM1_2 <- list(PBRM1KO_SMARCA4_2 = PBRM1KO_SMARCA4_2, WT_PBRM1_2 = WT_PBRM1_2)
overlap_PBRM1KO_SMARCA4_and_PBRM1_2_overlap <- overlap(Venn(overlap_PBRM1KO_SMARCA4_and_PBRM1_2))
overlap_PBRM1KO_SMARCA4_and_SMARCA4_2 <- list(WT_SMARCA4_2 = WT_SMARCA4_2, PBRM1KO_SMARCA4_2 = PBRM1KO_SMARCA4_2)
overlap_PBRM1KO_SMARCA4_and_SMARCA4_2_overlap <- overlap(Venn(overlap_PBRM1KO_SMARCA4_and_SMARCA4_2))

SMARCA4_and_PBRM1_not_PBRM1KO_SMARCA4_2 <- setdiff(overlap_SMARCA4_and_PBRM1_2_overlap, c(PBRM1KO_SMARCA4_2))


PBRM1KO_SMARCA4_and_PBRM1_not_SMARCA4_2 <- setdiff(overlap_PBRM1KO_SMARCA4_and_PBRM1_2_overlap, c(WT_SMARCA4_2))


PBRM1KO_SMARCA4_and_SMARCA4_not_PBRM1_2 <- setdiff(overlap_PBRM1KO_SMARCA4_and_SMARCA4_2_overlap, c(WT_PBRM1_2))



## final groups for publication

# PBRM1 non-specific kmers
PBRM1_non_specific_kmers <- c(overlap_SMARCA4_and_PBRM1_and_PBRM1KO_SMARCA4_2_overlap, PBRM1KO_SMARCA4_and_PBRM1_not_SMARCA4_2)
write.table(as.data.frame(PBRM1_non_specific_kmers), file="PBRM1_non-specific_kmers_fc2.txt", sep=" ", row.names=F, col.names = F, quote = F)

# PBRM1-specific kmers
PBRM1_specific_kmers <- c(Just_WT_PBRM1_2, SMARCA4_and_PBRM1_not_PBRM1KO_SMARCA4_2)
write.table(as.data.frame(PBRM1_specific_kmers), file="PBRM1_specific_kmers_fc2.txt", sep=" ", row.names=F, col.names = F, quote = F)

# KO-specific kmers

write.table(as.data.frame(Just_PBRM1KO_SMARCA4), file="KO_specific_kmers_fc2.txt", sep=" ", row.names=F, col.names = F, quote = F)
