setwd(system("pwd", intern = T) )


if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos = "http://cran.us.r-project.org")



install.packages("devEMF",version='4.5',repos = "http://cran.us.r-project.org")

require(GenomicRanges)
library(EnrichedHeatmap)
library(rtracklayer)
library(circlize)
library(data.table)
library(devEMF)

# prepare SMARCA4 heatmap from figure 3E

# make targets file from SMARCA4 peaks - centre around peak - 5kb either side

targets_Sm_cen <- makeGRangesFromDataFrame(
  df = fread("unique_reads/MACS2_peaks/overlap_2/CEN_peaks/WT_SMARCA4_vs_IgG_overlap_cen_peaks_merge.bed", header = FALSE, data.table = FALSE),
  seqnames.field = "V1", start.field = "V2", end.field = "V3")

ExtendSize <- 5000
targets_Sm_cen.extended  <- resize(targets_Sm_cen, fix = "center", width = ExtendSize*2)

# load bigwig files

# CUT&RUN

# WT SMARCA4

BigWig_Sm_3 <- rtracklayer::import("unique_reads/bigwig/WT_SMARCA4_LH_r3.bs10.bw", 
                                   format = "BigWig", 
                                   selection = BigWigSelection(targets_Sm_cen.extended))

normMatrix_Sm_3 <- normalizeToMatrix(signal = BigWig_Sm_3, 
                                     target = resize(targets_Sm_cen, fix = "center", width = 1), 
                                     background = 0, 
                                     keep = c(0, 0.99),      #/ minimal value to the 99th percentile
                                     target_ratio = 0,
                                     mean_mode = "coverage",       #/ see ?EnrichedHeatmap on other options
                                     value_column = "score", #/ = the name of the 4th column of the bigwig
                                     extend = ExtendSize)

col_fun_Sm = circlize::colorRamp2(quantile(normMatrix_Sm_3, c(0, .99)), c("white", "#00AD8A"))

# WT IgG

BigWig_Ig_3 <- rtracklayer::import("unique_reads/bigwig/WT_IgG_LrmDupH_r3.bs10.bw", 
                                   format = "BigWig", 
                                   selection = BigWigSelection(targets_Sm_cen.extended))

normMatrix_Ig_3 <- normalizeToMatrix(signal = BigWig_Ig_3, 
                                     target = resize(targets_Sm_cen, fix = "center", width = 1), 
                                     background = 0, 
                                     keep = c(0, 0.99),      #/ minimal value to the 99th percentile
                                     target_ratio = 0,
                                     mean_mode = "coverage",       #/ see ?EnrichedHeatmap on other options
                                     value_column = "score", #/ = the name of the 4th column of the bigwig
                                     extend = ExtendSize)

col_fun_Ig = circlize::colorRamp2(quantile(normMatrix_Ig_3, c(0, .99)), c("white", "#919191"))


# ChIP-seq

# GSE71510 - HCT116

# Parental SMARCA4

BigWig_P4 <- rtracklayer::import("unique_reads/ChIP_seq/bigwig/HCT116_Parental_SMARCA4.bs10.bw", 
                                 format = "BigWig", 
                                 selection = BigWigSelection(targets_Sm_cen.extended))

normMatrix_P4 <- normalizeToMatrix(signal = BigWig_P4, 
                                   target = resize(targets_Sm_cen, fix = "center", width = 1), 
                                   background = 0, 
                                   keep = c(0, 0.99),      #/ minimal value to the 99th percentile
                                   target_ratio = 0,
                                   mean_mode = "coverage",       #/ see ?EnrichedHeatmap on other options
                                   value_column = "score", #/ = the name of the 4th column of the bigwig
                                   extend = ExtendSize)

col_fun_P4 = circlize::colorRamp2(quantile(normMatrix_P4, c(0, .99)), c("white", "#00ab47"))

# Parental Input

BigWig_Pi <- rtracklayer::import("unique_reads/ChIP_seq/bigwig/HCT116_Parental_input.bs10.bw", 
                                 format = "BigWig", 
                                 selection = BigWigSelection(targets_Sm_cen.extended))

normMatrix_Pi <- normalizeToMatrix(signal = BigWig_Pi, 
                                   target = resize(targets_Sm_cen, fix = "center", width = 1), 
                                   background = 0, 
                                   keep = c(0, 0.99),      #/ minimal value to the 99th percentile
                                   target_ratio = 0,
                                   mean_mode = "coverage",       #/ see ?EnrichedHeatmap on other options
                                   value_column = "score", #/ = the name of the 4th column of the bigwig
                                   extend = ExtendSize)

col_fun_Pi = circlize::colorRamp2(quantile(normMatrix_Pi, c(0, .99)), c("white", "#919191"))



# GSE152681 - 786O

# Parental SMARCA4

BigWig_Sm_7 <- rtracklayer::import("unique_reads/ChIP_seq/bigwig/SMARCA4_786O_PBRM1WT.bs10.bw", 
                                   format = "BigWig", 
                                   selection = BigWigSelection(targets_Sm_cen.extended))

normMatrix_Sm_7 <- normalizeToMatrix(signal = BigWig_Sm_7, 
                                     target = resize(targets_Sm_cen, fix = "center", width = 1), 
                                     background = 0, 
                                     keep = c(0, 0.99),      #/ minimal value to the 99th percentile
                                     target_ratio = 0,
                                     mean_mode = "coverage",       #/ see ?EnrichedHeatmap on other options
                                     value_column = "score", #/ = the name of the 4th column of the bigwig
                                     extend = ExtendSize)

col_fun_Sm7 = circlize::colorRamp2(quantile(normMatrix_Sm_7, c(0, .99)), c("white", "#50A315"))

# Parental Input

BigWig_7i <- rtracklayer::import("unique_reads/ChIP_seq/bigwig/Input_SE_786O_PBRM1WT.bs10.bw", 
                                 format = "BigWig", 
                                 selection = BigWigSelection(targets_Sm_cen.extended))

normMatrix_7i <- normalizeToMatrix(signal = BigWig_7i, 
                                   target = resize(targets_Sm_cen, fix = "center", width = 1), 
                                   background = 0, 
                                   keep = c(0, 0.99),      #/ minimal value to the 99th percentile
                                   target_ratio = 0,
                                   mean_mode = "coverage",       #/ see ?EnrichedHeatmap on other options
                                   value_column = "score", #/ = the name of the 4th column of the bigwig
                                   extend = ExtendSize)

col_fun_7i = circlize::colorRamp2(quantile(normMatrix_7i, c(0, .99)), c("white", "#919191"))

# Prepare heatmap

EH_Sm <- EnrichedHeatmap( mat = normMatrix_Sm_3, name = "SMARCA4", 
                          pos_line = FALSE, #/ no dashed lines around the start
                          border = FALSE,   #/ no box around heatmap
                          col = col_fun_Sm,    #/ color gradients from above
                          column_title = "SMARCA4 rep3", #/ column title 
                          column_title_gp = gpar(fontsize = 15, fontfamily = "sans"),
                          use_raster = TRUE, raster_quality = 10, raster_device = "png",
                          #/ turn off background colors
                          rect_gp = gpar(col = "transparent"),  axis_name_gp = gpar(fontsize = 15), 
                          #/ legend options:
                          heatmap_legend_param = list(
                            labels_gp = gpar(font = 15),
                            title = "normalized counts"),
                          #/ options for the profile plot on top of the heatmap:
                          top_annotation = HeatmapAnnotation(
                            enriched = anno_enriched(
                              gp = gpar(col = 1, lty = 1, lwd=2), axis_param = list(gp = gpar(fontsize =15)),
                              col="black", ylim = c(0,30)))) +
  EnrichedHeatmap( mat = normMatrix_Ig_3, name = "IgG",
                   pos_line = FALSE, #/ no dashed lines around the start
                   border = FALSE,   #/ no box around heatmap
                   col = col_fun_Ig,    #/ color gradients from above
                   column_title = "IgG rep3", #/ column title 
                   column_title_gp = gpar(fontsize = 15, fontfamily = "sans"),
                   use_raster = TRUE, raster_quality = 10, raster_device = "png",
                   #/ turn off background colors
                   rect_gp = gpar(col = "transparent"),  axis_name_gp = gpar(fontsize = 15), 
                   #/ legend options:
                   heatmap_legend_param = list(
                     labels_gp = gpar(font = 2),
                     title = "normalized counts"),
                   #/ options for the profile plot on top of the heatmap:
                   top_annotation = HeatmapAnnotation(
                     enriched = anno_enriched(
                       gp = gpar(col = 1, lty = 1, lwd=2), axis_param = list(gp = gpar(fontsize =15)),
                       col="black", ylim = c(0,30)))) +
  EnrichedHeatmap( mat = normMatrix_P4, name = "HCT116 SMARCA4",
                   pos_line = FALSE, #/ no dashed lines around the start
                   border = FALSE,   #/ no box around heatmap
                   col = col_fun_P4,    #/ color gradients from above
                   column_title = "HCT116 parental SMARCA4", #/ column title 
                   column_title_gp = gpar(fontsize = 15, fontfamily = "sans"),
                   use_raster = TRUE, raster_quality = 10, raster_device = "png",
                   #/ turn off background colors
                   rect_gp = gpar(col = "transparent"),  axis_name_gp = gpar(fontsize = 15), 
                   #/ legend options:
                   heatmap_legend_param = list(
                     labels_gp = gpar(font = 15),
                     title = "normalized counts"),
                   #/ options for the profile plot on top of the heatmap:
                   top_annotation = HeatmapAnnotation(
                     enriched = anno_enriched(
                       gp = gpar(col = 1, lty = 1, lwd=2), axis_param = list(gp = gpar(fontsize =15)),
                       col="black", ylim = c(0,3)))) +
  EnrichedHeatmap( mat = normMatrix_Pi, name = "HCT116 input", 
                   pos_line = FALSE, #/ no dashed lines around the start
                   border = FALSE,   #/ no box around heatmap
                   col = col_fun_Pi,    #/ color gradients from above
                   column_title = "HCT116 Parental input", #/ column title 
                   column_title_gp = gpar(fontsize = 15, fontfamily = "sans"),
                   use_raster = TRUE, raster_quality = 10, raster_device = "png",
                   #/ turn off background colors
                   rect_gp = gpar(col = "transparent"),  axis_name_gp = gpar(fontsize = 15), 
                   #/ legend options:
                   heatmap_legend_param = list(
                     labels_gp = gpar(font = 15),
                     title = "normalized counts"),
                   #/ options for the profile plot on top of the heatmap:
                   top_annotation = HeatmapAnnotation(
                     enriched = anno_enriched(
                       gp = gpar(col = 1, lty = 1, lwd=2), axis_param = list(gp = gpar(fontsize =15)),
                       col="black", ylim = c(0,3))))  +
  EnrichedHeatmap( mat = normMatrix_Sm_7, name = "786O SMARCA4",
                   pos_line = FALSE, #/ no dashed lines around the start
                   border = FALSE,   #/ no box around heatmap
                   col = col_fun_Sm7,    #/ color gradients from above
                   column_title = "786O parental SMARCA4", #/ column title 
                   column_title_gp = gpar(fontsize = 15, fontfamily = "sans"),
                   use_raster = TRUE, raster_quality = 10, raster_device = "png",
                   #/ turn off background colors
                   rect_gp = gpar(col = "transparent"),  axis_name_gp = gpar(fontsize = 15), 
                   #/ legend options:
                   heatmap_legend_param = list(
                     labels_gp = gpar(font = 2),
                     title = "normalized counts"),
                   #/ options for the profile plot on top of the heatmap:
                   top_annotation = HeatmapAnnotation(
                     enriched = anno_enriched(
                       gp = gpar(col = 1, lty = 1, lwd=2), axis_param = list(gp = gpar(fontsize =15)),
                       col="black", ylim = c(0,100)))) +
  EnrichedHeatmap( mat = normMatrix_7i, name = "786O input", 
                   pos_line = FALSE, #/ no dashed lines around the start
                   border = FALSE,   #/ no box around heatmap
                   col = col_fun_7i,    #/ color gradients from above
                   column_title = "786O Parental input", #/ column title 
                   column_title_gp = gpar(fontsize = 15, fontfamily = "sans"),
                   use_raster = TRUE, raster_quality = 10, raster_device = "png",
                   #/ turn off background colors
                   rect_gp = gpar(col = "transparent"),  axis_name_gp = gpar(fontsize = 15), 
                   #/ legend options:
                   heatmap_legend_param = list(
                     labels_gp = gpar(font = 15),
                     title = "normalized counts"),
                   #/ options for the profile plot on top of the heatmap:
                   top_annotation = HeatmapAnnotation(
                     enriched = anno_enriched(
                       gp = gpar(col = 1, lty = 1, lwd=2), axis_param = list(gp = gpar(fontsize =15)),
                       col="black", ylim = c(0,100))))


#/ end of EnrichedHeatmap function

#/ Save as pdf or emf to disk:

pdf("unique_reads/Enriched_heatmaps/SMARCA4_cen_peaks_rep3_HCT116_786O_SMARCA4_clustered.pdf", width = 16, height = 10)

draw(EH_Sm,                                 #/ plot the heatmap from above 
     heatmap_legend_side = "right",     #/ we want the legend below the heatmap
     annotation_legend_side = "right",  #/ legend on the bottom side
     padding = unit(c(4, 4, 4, 4), "mm"), #/ some padding to avoid labels beyond plot borders
     ht_gap = unit(c(12, 12, 12, 12, 12), "mm")
)

dev.off()

emf("unique_reads/Enriched_heatmaps/SMARCA4_cen_peaks_rep3_HCT116_786O_SMARCA4_clustered.emf", units="in", width = 16, height = 10, coordDPI=300)

draw(EH_Sm,                                 #/ plot the heatmap from above 
     heatmap_legend_side = "right",     #/ we want the legend below the heatmap
     annotation_legend_side = "right",  #/ legend on the bottom side
     padding = unit(c(4, 4, 4, 4), "mm"), #/ some padding to avoid labels beyond plot borders
     ht_gap = unit(c(12, 12, 12, 12, 12), "mm")
)

dev.off()





# prepare SMARCA4 heatmap from figure 3E

# make targets file from SMARCA4 peaks - centre around peak - 5kb either side

targets_Pb_cen <- makeGRangesFromDataFrame(
  df = fread("unique_reads/MACS2_peaks/overlap_2/CEN_peaks/WT_PBRM1_vs_IgG_overlap_cen_peaks_merge.bed", header = FALSE, data.table = FALSE),
  seqnames.field = "V1", start.field = "V2", end.field = "V3")

ExtendSize <- 5000
targets_Pb_cen.extended  <- resize(targets_Pb_cen, fix = "center", width = ExtendSize*2)

# load bigwig files

# CUT&RUN

# WT PBRM1

BigWig_Pb_3 <- rtracklayer::import("unique_reads/bigwig/WT_PBRM1_LH_r3.bs10.bw",
                                   format = "BigWig",
                                   selection = BigWigSelection(targets_Pb_cen.extended))

normMatrix_Pb_3 <- normalizeToMatrix(signal = BigWig_Pb_3,
                                     target = resize(targets_Pb_cen, fix = "center", width = 1),
                                     background = 0,
                                     keep = c(0, 0.99),      #/ minimal value to the 99th percentile
                                     target_ratio = 0,
                                     mean_mode = "coverage",       #/ see ?EnrichedHeatmap on other options
                                     value_column = "score", #/ = the name of the 4th column of the bigwig
                                     extend = ExtendSize)

col_fun_Pb = circlize::colorRamp2(quantile(normMatrix_Pb_3, c(0, .99)), c("white", "#009ADE"))

# WT IgG

BigWig_Ig_3 <- rtracklayer::import("unique_reads/bigwig/WT_IgG_LrmDupH_r3.bs10.bw",
                                   format = "BigWig",
                                   selection = BigWigSelection(targets_Pb_cen.extended))

normMatrix_Ig_3 <- normalizeToMatrix(signal = BigWig_Ig_3,
                                     target = resize(targets_Pb_cen, fix = "center", width = 1),
                                     background = 0,
                                     keep = c(0, 0.99),      #/ minimal value to the 99th percentile
                                     target_ratio = 0,
                                     mean_mode = "coverage",       #/ see ?EnrichedHeatmap on other options
                                     value_column = "score", #/ = the name of the 4th column of the bigwig
                                     extend = ExtendSize)

col_fun_Ig = circlize::colorRamp2(quantile(normMatrix_Ig_3, c(0, .99)), c("white", "#919191"))


# ChIP-seq


# GSE152681 - 786O


# Parental PBRM1

BigWig_Pb_7 <- rtracklayer::import("unique_reads/ChIP_seq/bigwig/PBRM1_786O_PBRM1WT.bs10.bw",
                                   format = "BigWig",
                                   selection = BigWigSelection(targets_Pb_cen.extended))

normMatrix_Pb_7 <- normalizeToMatrix(signal = BigWig_Pb_7,
                                     target = resize(targets_Pb_cen, fix = "center", width = 1),
                                     background = 0,
                                     keep = c(0, 0.99),      #/ minimal value to the 99th percentile
                                     target_ratio = 0,
                                     mean_mode = "coverage",       #/ see ?EnrichedHeatmap on other options
                                     value_column = "score", #/ = the name of the 4th column of the bigwig
                                     extend = ExtendSize)

col_fun_Pb7 = circlize::colorRamp2(quantile(normMatrix_Pb_7, c(0, .99)), c("white", "#50A315"))

# Parental Input

BigWig_7i <- rtracklayer::import("unique_reads/ChIP_seq/bigwig/Input_SE_786O_PBRM1WT.bs10.bw",
                                 format = "BigWig",
                                 selection = BigWigSelection(targets_Pb_cen.extended))

normMatrix_7i <- normalizeToMatrix(signal = BigWig_7i,
                                   target = resize(targets_Pb_cen, fix = "center", width = 1),
                                   background = 0,
                                   keep = c(0, 0.99),      #/ minimal value to the 99th percentile
                                   target_ratio = 0,
                                   mean_mode = "coverage",       #/ see ?EnrichedHeatmap on other options
                                   value_column = "score", #/ = the name of the 4th column of the bigwig
                                   extend = ExtendSize)

col_fun_7i = circlize::colorRamp2(quantile(normMatrix_7i, c(0, .99)), c("white", "#919191"))

# GSE152681 - HK2

# Parental PBRM1

BigWig_Pb_HK2 <- rtracklayer::import("unique_reads/ChIP_seq/bigwig/PBRM1_HK2.bs10.bw", 
                                     format = "BigWig", 
                                     selection = BigWigSelection(targets_Pb_cen.extended))

normMatrix_Pb_HK2 <- normalizeToMatrix(signal = BigWig_Pb_HK2, 
                                       target = resize(targets_Pb_cen, fix = "center", width = 1), 
                                       background = 0, 
                                       keep = c(0, 0.99),      #/ minimal value to the 99th percentile
                                       target_ratio = 0,
                                       mean_mode = "coverage",       #/ see ?EnrichedHeatmap on other options
                                       value_column = "score", #/ = the name of the 4th column of the bigwig
                                       extend = ExtendSize)

col_fun_PbH = circlize::colorRamp2(quantile(normMatrix_Pb_HK2, c(0, .99)), c("white", "#3596E1"))

# Parental Input

BigWig_Hi <- rtracklayer::import("unique_reads/ChIP_seq/bigwig/Input_HK2.bs10.bw", 
                                 format = "BigWig", 
                                 selection = BigWigSelection(targets_Pb_cen.extended))

normMatrix_Hi <- normalizeToMatrix(signal = BigWig_Hi, 
                                   target = resize(targets_Pb_cen, fix = "center", width = 1), 
                                   background = 0, 
                                   keep = c(0, 0.99),      #/ minimal value to the 99th percentile
                                   target_ratio = 0,
                                   mean_mode = "coverage",       #/ see ?EnrichedHeatmap on other options
                                   value_column = "score", #/ = the name of the 4th column of the bigwig
                                   extend = ExtendSize)

col_fun_Hi = circlize::colorRamp2(quantile(normMatrix_Hi, c(0, .99)), c("white", "#919191"))


# GSE152681 - A498

# Parental PBRM1

BigWig_Pb_A <- rtracklayer::import("unique_reads/ChIP_seq/bigwig/PBRM1_A498.bs10.bw", 
                                   format = "BigWig", 
                                   selection = BigWigSelection(targets_Pb_cen.extended))

normMatrix_Pb_A <- normalizeToMatrix(signal = BigWig_Pb_A, 
                                     target = resize(targets_Pb_cen, fix = "center", width = 1), 
                                     background = 0, 
                                     keep = c(0, 0.99),      #/ minimal value to the 99th percentile
                                     target_ratio = 0,
                                     mean_mode = "coverage",       #/ see ?EnrichedHeatmap on other options
                                     value_column = "score", #/ = the name of the 4th column of the bigwig
                                     extend = ExtendSize)

col_fun_PbA = circlize::colorRamp2(quantile(normMatrix_Pb_A, c(0, .99)), c("white", "#00A4CE"))


# Parental Input

BigWig_Ai <- rtracklayer::import("unique_reads/ChIP_seq/bigwig/Input_A498.bs10.bw", 
                                 format = "BigWig", 
                                 selection = BigWigSelection(targets_Pb_cen.extended))

normMatrix_Ai <- normalizeToMatrix(signal = BigWig_Ai, 
                                   target = resize(targets_Pb_cen, fix = "center", width = 1), 
                                   background = 0, 
                                   keep = c(0, 0.99),      #/ minimal value to the 99th percentile
                                   target_ratio = 0,
                                   mean_mode = "coverage",       #/ see ?EnrichedHeatmap on other options
                                   value_column = "score", #/ = the name of the 4th column of the bigwig
                                   extend = ExtendSize)

col_fun_Ai = circlize::colorRamp2(quantile(normMatrix_Ai, c(0, .99)), c("white", "#919191"))


# Prepare heatmap


EH_Pb <- EnrichedHeatmap( mat = normMatrix_Pb_3, name = "PBRM1", 
                          pos_line = FALSE, #/ no dashed lines around the start
                          border = FALSE,   #/ no box around heatmap
                          col = col_fun_Pb,    #/ color gradients from above
                          column_title = "PBRM1 rep3", #/ column title 
                          column_title_gp = gpar(fontsize = 15, fontfamily = "sans"),
                          use_raster = TRUE, raster_quality = 10, raster_device = "png",
                          #/ turn off background colors
                          rect_gp = gpar(col = "transparent"),  axis_name_gp = gpar(fontsize = 15),
                          #/ legend options:
                          heatmap_legend_param = list(
                            labels_gp = gpar(font = 15),
                            title = "normalized counts"),
                          #/ options for the profile plot on top of the heatmap:
                          top_annotation = HeatmapAnnotation(
                            enriched = anno_enriched(
                              gp = gpar(col = 1, lty = 1, lwd=2), axis_param = list(gp = gpar(fontsize =15)),
                              col="black", ylim = c(0,15)))) +
  EnrichedHeatmap( mat = normMatrix_Ig_3, name = "IgG",
                   pos_line = FALSE, #/ no dashed lines around the start
                   border = FALSE,   #/ no box around heatmap
                   col = col_fun_Ig,    #/ color gradients from above
                   column_title = "IgG rep3", #/ column title 
                   column_title_gp = gpar(fontsize = 15, fontfamily = "sans"),
                   use_raster = TRUE, raster_quality = 10, raster_device = "png",
                   #/ turn off background colors
                   rect_gp = gpar(col = "transparent"),  axis_name_gp = gpar(fontsize = 15),
                   #/ legend options:
                   heatmap_legend_param = list(
                     labels_gp = gpar(font = 15),
                     title = "normalized counts"),
                   #/ options for the profile plot on top of the heatmap:
                   top_annotation = HeatmapAnnotation(
                     enriched = anno_enriched(
                       gp = gpar(col = 1, lty = 1, lwd=2), axis_param = list(gp = gpar(fontsize =15)),
                       col="black", ylim = c(0,15)))) +
  EnrichedHeatmap( mat = normMatrix_Pb_7, name = "786O PBRM1",
                   pos_line = FALSE, #/ no dashed lines around the start
                   border = FALSE,   #/ no box around heatmap
                   col = col_fun_Pb7,    #/ color gradients from above
                   column_title = "786O parental PBRM1", #/ column title 
                   column_title_gp = gpar(fontsize = 15, fontfamily = "sans"),
                   use_raster = TRUE, raster_quality = 10, raster_device = "png",
                   #/ turn off background colors
                   rect_gp = gpar(col = "transparent"),  axis_name_gp = gpar(fontsize = 15),
                   #/ legend options:
                   heatmap_legend_param = list(
                     labels_gp = gpar(font = 15),
                     title = "normalized counts"),
                   #/ options for the profile plot on top of the heatmap:
                   top_annotation = HeatmapAnnotation(
                     enriched = anno_enriched(
                       gp = gpar(col = 1, lty = 1, lwd=2), axis_param = list(gp = gpar(fontsize =15)),
                       col="black", ylim = c(0,40)))) +
  EnrichedHeatmap( mat = normMatrix_7i, name = "786O input", 
                   pos_line = FALSE, #/ no dashed lines around the start
                   border = FALSE,   #/ no box around heatmap
                   col = col_fun_7i,    #/ color gradients from above
                   column_title = "786O Parental input", #/ column title 
                   column_title_gp = gpar(fontsize = 15, fontfamily = "sans"),
                   use_raster = TRUE, raster_quality = 10, raster_device = "png",
                   #/ turn off background colors
                   rect_gp = gpar(col = "transparent"),  axis_name_gp = gpar(fontsize = 15),
                   #/ legend options:
                   heatmap_legend_param = list(
                     labels_gp = gpar(font = 15),
                     title = "normalized counts"),
                   #/ options for the profile plot on top of the heatmap:
                   top_annotation = HeatmapAnnotation(
                     enriched = anno_enriched(
                       gp = gpar(col = 1, lty = 1, lwd=2), axis_param = list(gp = gpar(fontsize =15)),
                       col="black", ylim = c(0,40))))  +
  EnrichedHeatmap( mat = normMatrix_Pb_HK2, name = "HK2 PBRM1",
                   pos_line = FALSE, #/ no dashed lines around the start
                   border = FALSE,   #/ no box around heatmap
                   col = col_fun_PbH,    #/ color gradients from above
                   column_title = "HK2 PBRM1", #/ column title 
                   column_title_gp = gpar(fontsize = 15, fontfamily = "sans"),
                   use_raster = TRUE, raster_quality = 10, raster_device = "png",
                   #/ turn off background colors
                   rect_gp = gpar(col = "transparent"), axis_name_gp = gpar(fontsize = 15), 
                   #/ legend options:
                   heatmap_legend_param = list(
                     labels_gp = gpar(font = 15),
                     title = "normalized counts"),
                   #/ options for the profile plot on top of the heatmap:
                   top_annotation = HeatmapAnnotation(
                     enriched = anno_enriched(
                       gp = gpar(col = 1, lty = 1, lwd=2), axis_param = list(gp = gpar(fontsize =15)),
                       col="black", ylim = c(0,20)))) +
  EnrichedHeatmap( mat = normMatrix_Hi, name = "HK2 input", 
                   pos_line = FALSE, #/ no dashed lines around the start
                   border = FALSE,   #/ no box around heatmap
                   col = col_fun_Hi,    #/ color gradients from above
                   column_title = "HK2 input", #/ column title 
                   column_title_gp = gpar(fontsize = 15, fontfamily = "sans"),
                   use_raster = TRUE, raster_quality = 10, raster_device = "png",
                   #/ turn off background colors
                   rect_gp = gpar(col = "transparent"),  axis_name_gp = gpar(fontsize = 15),
                   #/ legend options:
                   heatmap_legend_param = list(
                     labels_gp = gpar(font = 15),
                     title = "normalized counts"),
                   #/ options for the profile plot on top of the heatmap:
                   top_annotation = HeatmapAnnotation(
                     enriched = anno_enriched(
                       gp = gpar(col = 1, lty = 1, lwd=2), axis_param = list(gp = gpar(fontsize =15)),
                       col="black", ylim = c(0,20))))  +
  EnrichedHeatmap( mat = normMatrix_Pb_A, name = "A498 PBRM1",
                   pos_line = FALSE, #/ no dashed lines around the start
                   border = FALSE,   #/ no box around heatmap
                   col = col_fun_PbA,    #/ color gradients from above
                   column_title = "A498 PBRM1", #/ column title 
                   column_title_gp = gpar(fontsize = 15, fontfamily = "sans"),
                   use_raster = TRUE, raster_quality = 10, raster_device = "png",
                   #/ turn off background colors
                   rect_gp = gpar(col = "transparent"),  axis_name_gp = gpar(fontsize = 15),
                   #/ legend options:
                   heatmap_legend_param = list(
                     labels_gp = gpar(font = 15),
                     title = "normalized counts"),
                   #/ options for the profile plot on top of the heatmap:
                   top_annotation = HeatmapAnnotation(
                     enriched = anno_enriched(
                       gp = gpar(col = 1, lty = 1, lwd=2), axis_param = list(gp = gpar(fontsize =15)),
                       col="black", ylim = c(0,15)))) +
  EnrichedHeatmap( mat = normMatrix_Ai, name = "A498 input", 
                   pos_line = FALSE, #/ no dashed lines around the start
                   border = FALSE,   #/ no box around heatmap
                   col = col_fun_Ai,    #/ color gradients from above
                   column_title = "A498 input", #/ column title 
                   column_title_gp = gpar(fontsize = 15, fontfamily = "sans"),
                   use_raster = TRUE, raster_quality = 10, raster_device = "png",
                   #/ turn off background colors
                   rect_gp = gpar(col = "transparent"),  axis_name_gp = gpar(fontsize = 15),
                   #/ legend options:
                   heatmap_legend_param = list(
                     labels_gp = gpar(font = 15),
                     title = "normalized counts"),
                   #/ options for the profile plot on top of the heatmap:
                   top_annotation = HeatmapAnnotation(
                     enriched = anno_enriched(
                       gp = gpar(col = 1, lty = 1, lwd=2), axis_param = list(gp = gpar(fontsize =15)),
                       col="black", ylim = c(0,15))))


#/ end of EnrichedHeatmap function

#/ Save as pdf to disk:
pdf("unique_reads/Enriched_heatmaps/PBRM1_cen_peaks_rep3_786O_HK2_A498_clustered.pdf", width = 18, height = 10)

draw(EH_Pb,                                 #/ plot the heatmap from above 
     heatmap_legend_side = "right",     #/ we want the legend below the heatmap
     annotation_legend_side = "right",  #/ legend on the bottom side
     padding = unit(c(4, 4, 4, 4), "mm"), #/ some padding to avoid labels beyond plot borders
     ht_gap = unit(c(12, 12, 12, 12, 12, 12, 12), "mm")
)

dev.off()

emf("unique_reads/Enriched_heatmaps/PBRM1_cen_peaks_rep3_786O_HK2_A498_clustered.emf", units="in", width = 18, height = 10, coordDPI=300)

draw(EH_Pb,                                 #/ plot the heatmap from above 
     heatmap_legend_side = "right",     #/ we want the legend below the heatmap
     annotation_legend_side = "right",  #/ legend on the bottom side
     padding = unit(c(4, 4, 4, 4), "mm"), #/ some padding to avoid labels beyond plot borders
     ht_gap = unit(c(12, 12, 12, 12, 12, 12, 12), "mm")
)

dev.off()



