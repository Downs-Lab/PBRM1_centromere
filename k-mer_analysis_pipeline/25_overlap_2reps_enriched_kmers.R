
#set file path to current directory

setwd(system("pwd", intern = T) )

install.packages("ggVennDiagram",version='1.5.2',repos = "http://cran.us.r-project.org")
install.packages("RVenn",version='1.1.0',repos = "http://cran.us.r-project.org")


library(ggVennDiagram)
library(ggplot2)
library(RVenn)
library(R.utils)

args <- commandArgs(TRUE)


Sample_name = args[1]


Sample_1_2 <- read.delim(file =(paste("kmer_centromere/enriched_kmers/",Sample_name,"_1_vs_Control_fc2.fa.txt", 
                                          sep="")), header = F, sep = "")
Sample_1_2 <- Sample_1_2[,1]

Sample_2_2 <- read.delim(file =(paste("kmer_centromere/enriched_kmers/",Sample_name,"_2_vs_Control_fc2.fa.txt", 
                                      sep="")), header = F, sep = "")
Sample_2_2 <- Sample_2_2[,1]

Sample_3_2 <- read.delim(file =(paste("kmer_centromere/enriched_kmers/",Sample_name,"_3_vs_Control_fc2.fa.txt", 
                                      sep="")), header = F, sep = "")
Sample_3_2 <- Sample_3_2[,1]



Sample_1_name <- (paste(Sample_name,"_1", sep=""))
Sample_2_name <- (paste(Sample_name,"_2", sep=""))
Sample_3_name <- (paste(Sample_name,"_3", sep=""))

venn_Sample_list <- list(Sample_1_2 = Sample_1_2, Sample_2_2 = Sample_2_2, Sample_3_2 = Sample_3_2)

names(venn_Sample_list) <- c(Sample_1_name, Sample_2_name, Sample_3_name)

pdf(file=paste("kmer_centromere/enriched_kmers/overlap_2/",Sample_name,"_reps_venn.pdf", sep=""), width=5, height=5)
ggVennDiagram(venn_Sample_list, label = "count", label_alpha = 0) + scale_fill_gradient(low="white",high = "red") + 
  scale_color_manual(values = c(Sample_1_name = "black", Sample_2_name = "black", Sample_3_name ="black"))
dev.off()


# extract out enriched kmers that are found in at least 2 reps


Sample_1_2_2 <- overlap(Venn(venn_Sample_list), slice = c(Sample_1_name, Sample_2_name))
Sample_1_3_2 <- overlap(Venn(venn_Sample_list), slice = c(Sample_1_name, Sample_3_name))
Sample_2_3_2 <- overlap(Venn(venn_Sample_list), slice = c(Sample_2_name, Sample_3_name))
Sample_2_overlap <- append(Sample_1_2_2, Sample_1_3_2)
Sample_2_overlap <- append(Sample_2_overlap, Sample_2_3_2)
Sample_2_overlap <- paste(unique(Sample_2_overlap))
write.table(as.data.frame(Sample_2_overlap), file=paste("kmer_centromere/enriched_kmers/overlap_2/",Sample_name,"_overlap_2.txt", sep=""), sep=" ", row.names=F, col.names = F, quote = F)






