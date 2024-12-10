setwd(system("pwd", intern = T) )


if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos = "http://cran.us.r-project.org")



install.packages("devEMF",version='4.5',repos = "http://cran.us.r-project.org")


library(rtracklayer)
library(GenomeInfoDb)
library(AnnotationDbi)
library(ChIPseeker)
library(ensembldb)
library(org.Hs.eg.db)
library(devEMF)


# now you can reload the database without redoing earlier steps from 33_TxDb

txdb <- loadDb("genomes/txdb_t2t")


## Example ChIPseeker script below used to make a figure in the publication:

files_1 <- c("unique_reads/MACS2_peaks/overlap_2/WT_PBRM1_vs_IgG_overlap_peaks_merge.bed", "unique_reads/MACS2_peaks/overlap_2/WT_SMARCA4_vs_IgG_overlap_peaks_merge.bed", "unique_reads/MACS2_peaks/overlap_2/PBRM1KO_SMARCA4_vs_IgG_overlap_peaks_merge.bed")

names(files_1) <- c("PBRM1", "SMARCA4", "PBRM1-KO SMARCA4")


peakAnnoList <- lapply(files_1, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), annoDb="org.Hs.eg.db", verbose=FALSE)

# extract files as tiff

tiff('unique_reads/ChIPseeker/PBRM1_SMARCA4_PKO_peaks_merge_AnnoBar.tiff', units="in", width=8, height=4, res=300, compression = 'lzw')
plotAnnoBar(peakAnnoList)
dev.off()

tiff('unique_reads/ChIPseeker/PBRM1_SMARCA4_PKO_peaks_merge_PlotDistToTSS.tiff', units="in", width=8, height=4, res=300, compression = 'lzw')
plotDistToTSS(peakAnnoList)
dev.off()

# extract files in emf format

emf('unique_reads/ChIPseeker/PBRM1_SMARCA4_PKO_peaks_merge_AnnoBar.emf', units="in", width=8, height=4, coordDPI=300)
plotAnnoBar(peakAnnoList)
dev.off()

emf('unique_reads/ChIPseeker/PBRM1_SMARCA4_PKO_peaks_merge_PlotDistToTSS.emf', units="in", width=8, height=4, coordDPI=300)
plotDistToTSS(peakAnnoList)
dev.off()


# Equivalent example script with CEN and peri-CEN peaks


cen_files_1 <- c("unique_reads/MACS2_peaks/overlap_2/CEN_peaks/WT_PBRM1_vs_IgG_overlap_peaks_merge.bed", "unique_reads/MACS2_peaks/overlap_2/CEN_peaks/WT_SMARCA4_vs_IgG_overlap_peaks_merge.bed", "unique_reads/MACS2_peaks/overlap_2/CEN_peaks/PBRM1KO_SMARCA4_vs_IgG_overlap_peaks_merge.bed")

names(cen_files_1) <- c("PBRM1", "SMARCA4", "PBRM1-KO SMARCA4")


peakAnnoList_cen <- lapply(cen_files_1, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), annoDb="org.Hs.eg.db", verbose=FALSE)

# extract files as tiff

tiff('unique_reads/ChIPseeker/PBRM1_SMARCA4_PKO_cen_peaks_merge_AnnoBar.tiff', units="in", width=8, height=4, res=300, compression = 'lzw')
plotAnnoBar(peakAnnoList_cen)
dev.off()

tiff('unique_reads/ChIPseeker/PBRM1_SMARCA4_PKO_cen_peaks_merge_PlotDistToTSS.tiff', units="in", width=8, height=4, res=300, compression = 'lzw')
plotDistToTSS(peakAnnoList_cen)
dev.off()

# extract files in emf format

emf('unique_reads/ChIPseeker/PBRM1_SMARCA4_PKO_cen_peaks_merge_AnnoBar.emf', units="in", width=8, height=4, coordDPI=300)
plotAnnoBar(peakAnnoList_cen)
dev.off()

emf('unique_reads/ChIPseeker/PBRM1_SMARCA4_PKO_cen_peaks_merge_PlotDistToTSS.emf', units="in", width=8, height=4, coordDPI=300)
plotDistToTSS(peakAnnoList_cen)
dev.off()
