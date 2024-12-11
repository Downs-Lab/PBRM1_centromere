setwd(system("pwd", intern = T) )


install.packages("ggVennDiagram",version='1.5.2',repos = "http://cran.us.r-project.org")
install.packages("RVenn",version='1.1.0',repos = "http://cran.us.r-project.org")


library(ggVennDiagram)
library(ggplot2)
library(RVenn)
library(R.utils)


#load k-mer files

SMARCA4 <- read.delim(file ="kmer_centromere/enriched_kmers/overlap_2/WT_SMARCA4_overlap_2.txt",
                                          header = F, sep = "")
SMARCA4 <- SMARCA4[,1]


PBRM1 <- read.delim(file ="kmer_centromere/enriched_kmers/overlap_2/WT_PBRM1_overlap_2.txt",
                                          header = F, sep = "")
PBRM1 <- PBRM1[,1]


# prepare venn diagram

venn_Sample_list <- list(SMARCA4 = SMARCA4, PBRM1 = PBRM1)

pdf(file="kmer_centromere/enriched_kmers/overlap_2/venns/SMARCA4_PBRM1_venn.pdf", width=5, height=5)
ggVennDiagram(venn_Sample_list, label = "count", label_alpha = 0) + scale_fill_gradient(low="white",high = "red") + 
  scale_color_manual(values = c(SMARCA4 = "black", PBRM1 = "black"))
dev.off()
