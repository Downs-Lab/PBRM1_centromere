setwd(system("pwd", intern = T) )

install.packages("ggVennDiagram",version='1.5.2',repos = "http://cran.us.r-project.org")
install.packages("RVenn",version='1.1.0',repos = "http://cran.us.r-project.org")
install.packages("devEMF",version='4.5',repos = "http://cran.us.r-project.org")

library(ggVennDiagram)
library(ggplot2)
library(RVenn)
library(devEMF)

## Example script from publication data for parental and PBRM1-KO RPE1 H3K9me2 peak intersection
## This was the script used to create the final venn diagrams using the numbers generated with the intervene script


# WT_H3K9me2 and PBRM1KO_H3K9me2 peaks

# centromere peaks 

venn_WT_H3K9me2_and_PBRM1KO_H3K9me2 <- list(WT_H3K9me2 = 1:1923, PBRM1KO_H3K9me2 = c(1:268, 2001:2421))
  
ggVenn_WT_H3K9me2_and_PBRM1KO_H3K9me2 <- ggVennDiagram(venn_WT_H3K9me2_and_PBRM1KO_H3K9me2, label = "count", 
                                            label_alpha = 0, category.names = c("WT me2", 
                                                                                "KO me2"), set_size = 4) + 
scale_fill_gradient(low="white",high = "darkorange2") + 
scale_color_manual(values = c(WT_H3K9me2 = "black",PBRM1KO_H3K9me2 = "black")) + 
scale_x_continuous(expand = expansion(mult = .2))
  
  
emf('unique_reads/custom_venn/WT_H3K9me2_and_PBRM1KO_H3K9me2_CEN_peaks_venn.emf', units="in", width=6, height=6, coordDPI=300)
print(ggVenn_WT_H3K9me2_and_PBRM1KO_H3K9me2)
dev.off()

ggsave(file="unique_reads/custom_venn/WT_H3K9me2_and_PBRM1KO_H3K9me2_CEN_peaks_venn.svg", plot=ggVenn_WT_H3K9me2_and_PBRM1KO_H3K9me2, width=10, height=8)

# whole genome peaks

venn_WT_H3K9me2_and_PBRM1KO_H3K9me2 <- list(WT_H3K9me2 = 1:70912, PBRM1KO_H3K9me2 = c(1:12294, 80001:94338))

ggVenn_WT_H3K9me2_and_PBRM1KO_H3K9me2 <- ggVennDiagram(venn_WT_H3K9me2_and_PBRM1KO_H3K9me2, label = "count", 
                                          label_alpha = 0, category.names = c("WT me2", 
                                                                              "KO me2"), set_size = 4) + 
  scale_fill_gradient(low="white",high = "darkorange2") + 
  scale_color_manual(values = c(WT_H3K9me2 = "black",PBRM1KO_H3K9me2 = "black")) + 
  scale_x_continuous(expand = expansion(mult = .2))


emf('unique_reads/custom_venn/WT_H3K9me2_and_PBRM1KO_H3K9me2_genome_peaks_venn.emf', units="in", width=6, height=6, coordDPI=300)
print(ggVenn_WT_H3K9me2_and_PBRM1KO_H3K9me2)
dev.off()

ggsave(file="unique_reads/custom_venn/WT_H3K9me2_and_PBRM1KO_H3K9me2_genome_peaks_venn.svg", plot=ggVenn_WT_H3K9me2_and_PBRM1KO_H3K9me2, width=10, height=8)
