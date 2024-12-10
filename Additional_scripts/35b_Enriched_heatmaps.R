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

# prepare heatmaps from figure 4D

# PBRM1 specific peaks

# make targets file from PBRM1 specific peaks - centre around peak - 5kb either side

targets_Pb_spec <- makeGRangesFromDataFrame(
  df = fread("unique_reads/MACS2_peaks/overlap_2/CEN_peaks/PBRM1_specific_cen_peaks_merge.bed", header = FALSE, data.table = FALSE),
  seqnames.field = "V1", start.field = "V2", end.field = "V3")

ExtendSize <- 5000
targets_Pb_spec.extended  <- resize(targets_Pb_spec, fix = "center", width = ExtendSize*2)


# load bigwig files

# CUT&RUN

# WT SMARCA4

BigWig_Sm_3 <- rtracklayer::import("unique_reads/bigwig/WT_SMARCA4_LH_r3.bs10.bw", 
                                   format = "BigWig", 
                                   selection = BigWigSelection(targets_Pb_spec.extended))

normMatrix_Sm_3 <- normalizeToMatrix(signal = BigWig_Sm_3, 
                                     target = resize(targets_Pb_spec, fix = "center", width = 1), 
                                     background = 0, 
                                     keep = c(0, 0.99),      #/ minimal value to the 99th percentile
                                     target_ratio = 0,
                                     mean_mode = "coverage",       #/ see ?EnrichedHeatmap on other options
                                     value_column = "score", #/ = the name of the 4th column of the bigwig
                                     extend = ExtendSize)

col_fun_Sm = circlize::colorRamp2(quantile(normMatrix_Sm_3, c(0, .99)), c("white", "#00AD8A"))

# WT PBRM1

BigWig_Pb_3 <- rtracklayer::import("unique_reads/bigwig/WT_PBRM1_LH_r3.bs10.bw",
                                   format = "BigWig",
                                   selection = BigWigSelection(targets_Pb_spec.extended))

normMatrix_Pb_3 <- normalizeToMatrix(signal = BigWig_Pb_3,
                                     target = resize(targets_Pb_spec, fix = "center", width = 1),
                                     background = 0,
                                     keep = c(0, 0.99),      #/ minimal value to the 99th percentile
                                     target_ratio = 0,
                                     mean_mode = "coverage",       #/ see ?EnrichedHeatmap on other options
                                     value_column = "score", #/ = the name of the 4th column of the bigwig
                                     extend = ExtendSize)

col_fun_Pb = circlize::colorRamp2(quantile(normMatrix_Pb_3, c(0, .99)), c("white", "#009ADE"))

#PBRM1-KO SMARCA4

BigWig_PKO_3 <- rtracklayer::import("unique_reads/bigwig/PBRM1KO_SMARCA4_LH_r3.bs10.bw", 
                                    format = "BigWig", 
                                    selection = BigWigSelection(targets_Pb_spec.extended))

normMatrix_PKO_3 <- normalizeToMatrix(signal = BigWig_PKO_3, 
                                      target = resize(targets_Pb_spec, fix = "center", width = 1), 
                                      background = 0, 
                                      keep = c(0, 0.99),      #/ minimal value to the 99th percentile
                                      target_ratio = 0,
                                      mean_mode = "coverage",       #/ see ?EnrichedHeatmap on other options
                                      value_column = "score", #/ = the name of the 4th column of the bigwig
                                      extend = ExtendSize)

col_fun_PKO = circlize::colorRamp2(quantile(normMatrix_PKO_3, c(0, .99)), c("white", "#50A315"))

# WT IgG

BigWig_Ig_3 <- rtracklayer::import("unique_reads/bigwig/WT_IgG_LrmDupH_r3.bs10.bw",
                                   format = "BigWig",
                                   selection = BigWigSelection(targets_Pb_spec.extended))

normMatrix_Ig_3 <- normalizeToMatrix(signal = BigWig_Ig_3,
                                     target = resize(targets_Pb_spec, fix = "center", width = 1),
                                     background = 0,
                                     keep = c(0, 0.99),      #/ minimal value to the 99th percentile
                                     target_ratio = 0,
                                     mean_mode = "coverage",       #/ see ?EnrichedHeatmap on other options
                                     value_column = "score", #/ = the name of the 4th column of the bigwig
                                     extend = ExtendSize)

col_fun_Ig = circlize::colorRamp2(quantile(normMatrix_Ig_3, c(0, .99)), c("white", "#919191"))


# Prepare heatmap

#### PBRM1 specific cen peaks - SWI/SNF signal 


set.seed(123)

EH_PS <- 
  EnrichedHeatmap(mat = normMatrix_Pb_3, name = "PBRM1",
                  pos_line = FALSE, #/ no dashed lines around the start
                  border = FALSE,   #/ no box around heatmap
                  col = col_fun_Pb,    #/ color gradients from above
                  column_title = "PBRM1 rep3", #/ column title 
                  column_title_gp = gpar(fontsize = 15, fontfamily = "sans"),
                  use_raster = TRUE, raster_quality = 10, raster_device = "png",
                  #/ turn off background colors
                  rect_gp = gpar(col = "transparent"), axis_name_gp = gpar(fontsize = 15),
                  #/ legend options:
                  heatmap_legend_param = list(
                    labels_gp = gpar(font = 15),
                    title = "normalized counts"),
                  #/ options for the profile plot on top of the heatmap:
                  top_annotation = HeatmapAnnotation(enriched = anno_enriched(
                    gp = gpar(col = 1:4, lty = 1:4, lwd=2), axis_param = list(gp = gpar(fontsize =15)),
                    col="black", ylim = c(0,10)))) +
  EnrichedHeatmap(mat = normMatrix_Sm_3, name = "SMARCA4", 
                  pos_line = FALSE, #/ no dashed lines around the start
                  border = FALSE,   #/ no box around heatmap
                  col = col_fun_Sm,    #/ color gradients from above
                  column_title = "SMARCA4 rep3", #/ column title 
                  column_title_gp = gpar(fontsize = 15, fontfamily = "sans"),
                  use_raster = TRUE, raster_quality = 10, raster_device = "png",
                  #/ turn off background colors
                  rect_gp = gpar(col = "transparent"), axis_name_gp = gpar(fontsize = 15),
                  #/ legend options:
                  heatmap_legend_param = list(
                    labels_gp = gpar(font = 15),
                    title = "normalized counts"),
                  #/ options for the profile plot on top of the heatmap:
                  top_annotation = HeatmapAnnotation(enriched = anno_enriched(
                    gp = gpar(col = 1:4, lty = 1:4, lwd=2), axis_param = list(gp = gpar(fontsize =15)),
                    col="black", ylim = c(0,10)))) +
  EnrichedHeatmap(mat = normMatrix_PKO_3, name = "PBRM1_KO SMARCA4",
                  pos_line = FALSE, #/ no dashed lines around the start
                  border = FALSE,   #/ no box around heatmap
                  col = col_fun_PKO,    #/ color gradients from above
                  column_title = "PBRM1-KO SMARCA4 rep3", #/ column title 
                  column_title_gp = gpar(fontsize = 15, fontfamily = "sans"),
                  use_raster = TRUE, raster_quality = 10, raster_device = "png",
                  #/ turn off background colors
                  rect_gp = gpar(col = "transparent"), axis_name_gp = gpar(fontsize = 15),
                  #/ legend options:
                  heatmap_legend_param = list(
                    labels_gp = gpar(font = 15),
                    title = "normalized counts"),
                  #/ options for the profile plot on top of the heatmap:
                  top_annotation = HeatmapAnnotation(enriched = anno_enriched(
                    gp = gpar(col = 1:4, lty = 1:4, lwd=2), axis_param = list(gp = gpar(fontsize =15)),
                    col="black", ylim = c(0,10)))) +
  EnrichedHeatmap(mat = normMatrix_Ig_3, name = "IgG",
                  pos_line = FALSE, #/ no dashed lines around the start
                  border = FALSE,   #/ no box around heatmap
                  col = col_fun_Ig,    #/ color gradients from above
                  column_title = "IgG rep3", #/ column title 
                  column_title_gp = gpar(fontsize = 15, fontfamily = "sans"),
                  use_raster = TRUE, raster_quality = 10, raster_device = "png",
                  #/ turn off background colors
                  rect_gp = gpar(col = "transparent"), axis_name_gp = gpar(fontsize = 15),
                  #/ legend options:
                  heatmap_legend_param = list(
                    labels_gp = gpar(font = 15),
                    title = "normalized counts"),
                  #/ options for the profile plot on top of the heatmap:
                  top_annotation = HeatmapAnnotation(enriched = anno_enriched(
                    gp = gpar(col = 1:4, lty = 1:4, lwd=2), axis_param = list(gp = gpar(fontsize =15)),
                    col="black", ylim = c(0,10))))


# save pdf, tiff and emf files to disk

pdf("unique_reads/Enriched_heatmaps/PBRM1_specific_cen_peaks_rep3_SWI_SNF.pdf", width = 15, height = 12)


draw(EH_PS,                                 #/ plot the heatmap from above 
     heatmap_legend_side = "right",     #/ we want the legend below the heatmap
     annotation_legend_side = "right",  #/ legend on the bottom side
     padding = unit(c(4, 4, 4, 4), "mm"), #/ some padding to avoid labels beyond plot borders
     ht_gap = unit(c(10, 10, 10), "mm"))



dev.off()



tiff('unique_reads/Enriched_heatmaps/PBRM1_specific_cen_peaks_rep3_SWI_SNF.tiff', units="in", width=12, height=12, res=300, compression = 'lzw')

draw(EH_PS,                                 #/ plot the heatmap from above 
     heatmap_legend_side = "right",     #/ we want the legend below the heatmap
     annotation_legend_side = "right",  #/ legend on the bottom side
     padding = unit(c(4, 4, 4, 4), "mm"), #/ some padding to avoid labels beyond plot borders
     ht_gap = unit(c(10, 10, 10), "mm"))



dev.off()



emf('unique_reads/Enriched_heatmaps/PBRM1_specific_cen_peaks_rep3_SWI_SNF.emf', units="in", width=6, height=6, coordDPI=300)

draw(EH_PS,                                 #/ plot the heatmap from above 
     heatmap_legend_side = "right",     #/ we want the legend below the heatmap
     annotation_legend_side = "right",  #/ legend on the bottom side
     padding = unit(c(4, 4, 4, 4), "mm"), #/ some padding to avoid labels beyond plot borders
     ht_gap = unit(c(10, 10, 10), "mm"))



dev.off()



# PBRM1 non-specific peaks

# make targets file from PBRM1 non-specific peaks - centre around peak - 5kb either side

targets_Pb_non <- makeGRangesFromDataFrame(
  df = fread("unique_reads/MACS2_peaks/overlap_2/CEN_peaks/PBRM1_specific_cen_peaks_merge.bed", header = FALSE, data.table = FALSE),
  seqnames.field = "V1", start.field = "V2", end.field = "V3")

ExtendSize <- 5000
targets_Pb_non.extended  <- resize(targets_Pb_non, fix = "center", width = ExtendSize*2)


# load bigwig files

# CUT&RUN

# WT SMARCA4

BigWig_Sm_3 <- rtracklayer::import("unique_reads/bigwig/WT_SMARCA4_LH_r3.bs10.bw",
                                   format = "BigWig",
                                   selection = BigWigSelection(targets_Pb_non.extended))

normMatrix_Sm_3 <- normalizeToMatrix(signal = BigWig_Sm_3,
                                     target = resize(targets_Pb_non, fix = "center", width = 1),
                                     background = 0,
                                     keep = c(0, 0.99),      #/ minimal value to the 99th percentile
                                     target_ratio = 0,
                                     mean_mode = "coverage",       #/ see ?EnrichedHeatmap on other options
                                     value_column = "score", #/ = the name of the 4th column of the bigwig
                                     extend = ExtendSize)

col_fun_Sm = circlize::colorRamp2(quantile(normMatrix_Sm_3, c(0, .99)), c("white", "#00AD8A"))

# WT PBRM1

BigWig_Pb_3 <- rtracklayer::import("unique_reads/bigwig/WT_PBRM1_LH_r3.bs10.bw",
                                   format = "BigWig",
                                   selection = BigWigSelection(targets_Pb_non.extended))

normMatrix_Pb_3 <- normalizeToMatrix(signal = BigWig_Pb_3,
                                     target = resize(targets_Pb_non, fix = "center", width = 1),
                                     background = 0,
                                     keep = c(0, 0.99),      #/ minimal value to the 99th percentile
                                     target_ratio = 0,
                                     mean_mode = "coverage",       #/ see ?EnrichedHeatmap on other options
                                     value_column = "score", #/ = the name of the 4th column of the bigwig
                                     extend = ExtendSize)

col_fun_Pb = circlize::colorRamp2(quantile(normMatrix_Pb_3, c(0, .99)), c("white", "#009ADE"))

#PBRM1-KO SMARCA4

BigWig_PKO_3 <- rtracklayer::import("unique_reads/bigwig/PBRM1KO_SMARCA4_LH_r3.bs10.bw",
                                    format = "BigWig",
                                    selection = BigWigSelection(targets_Pb_non.extended))

normMatrix_PKO_3 <- normalizeToMatrix(signal = BigWig_PKO_3,
                                      target = resize(targets_Pb_non, fix = "center", width = 1),
                                      background = 0,
                                      keep = c(0, 0.99),      #/ minimal value to the 99th percentile
                                      target_ratio = 0,
                                      mean_mode = "coverage",       #/ see ?EnrichedHeatmap on other options
                                      value_column = "score", #/ = the name of the 4th column of the bigwig
                                      extend = ExtendSize)

col_fun_PKO = circlize::colorRamp2(quantile(normMatrix_PKO_3, c(0, .99)), c("white", "#50A315"))

# WT IgG

BigWig_Ig_3 <- rtracklayer::import("unique_reads/bigwig/WT_IgG_LrmDupH_r3.bs10.bw",
                                   format = "BigWig",
                                   selection = BigWigSelection(targets_Pb_non.extended))

normMatrix_Ig_3 <- normalizeToMatrix(signal = BigWig_Ig_3,
                                     target = resize(targets_Pb_non, fix = "center", width = 1),
                                     background = 0,
                                     keep = c(0, 0.99),      #/ minimal value to the 99th percentile
                                     target_ratio = 0,
                                     mean_mode = "coverage",       #/ see ?EnrichedHeatmap on other options
                                     value_column = "score", #/ = the name of the 4th column of the bigwig
                                     extend = ExtendSize)

col_fun_Ig = circlize::colorRamp2(quantile(normMatrix_Ig_3, c(0, .99)), c("white", "#919191"))


# Prepare heatmap

#### PBRM1 non-specific cen peaks - SWI/SNF signal

set.seed(123)

EH_NS <- 
  EnrichedHeatmap(mat = normMatrix_Pb_3, name = "PBRM1",
                  pos_line = FALSE, #/ no dashed lines around the start
                  border = FALSE,   #/ no box around heatmap
                  col = col_fun_Pb,    #/ color gradients from above
                  column_title = "PBRM1 rep3", #/ column title 
                  column_title_gp = gpar(fontsize = 12, fontfamily = "sans"),
                  use_raster = TRUE, raster_quality = 10, raster_device = "png",
                  #/ turn off background colors
                  rect_gp = gpar(col = "transparent"), axis_name_gp = gpar(fontsize = 15),
                  #/ legend options:
                  heatmap_legend_param = list(
                    labels_gp = gpar(font = 15),
                    title = "normalized counts"),
                  #/ options for the profile plot on top of the heatmap:
                  top_annotation = HeatmapAnnotation(enriched = anno_enriched(
                    gp = gpar(col = 1:4, lty = 1:4, lwd=2), axis_param = list(gp = gpar(fontsize =15)),
                    col="black", ylim = c(0,20)))) +
  EnrichedHeatmap(mat = normMatrix_Sm_3, name = "SMARCA4", 
                  pos_line = FALSE, #/ no dashed lines around the start
                  border = FALSE,   #/ no box around heatmap
                  col = col_fun_Sm,    #/ color gradients from above
                  column_title = "SMARCA4 rep3", #/ column title 
                  column_title_gp = gpar(fontsize = 12, fontfamily = "sans"),
                  use_raster = TRUE, raster_quality = 10, raster_device = "png",
                  #/ turn off background colors
                  rect_gp = gpar(col = "transparent"), axis_name_gp = gpar(fontsize = 15),
                  #/ legend options:
                  heatmap_legend_param = list(
                    labels_gp = gpar(font = 15),
                    title = "normalized counts"),
                  #/ options for the profile plot on top of the heatmap:
                  top_annotation = HeatmapAnnotation(enriched = anno_enriched(
                    gp = gpar(col = 1:4, lty = 1:4, lwd=2), axis_param = list(gp = gpar(fontsize =15)),
                    col="black", ylim = c(0,20)))) +
  EnrichedHeatmap(mat = normMatrix_PKO_3, name = "PBRM1_KO SMARCA4",
                  pos_line = FALSE, #/ no dashed lines around the start
                  border = FALSE,   #/ no box around heatmap
                  col = col_fun_PKO,    #/ color gradients from above
                  column_title = "PBRM1-KO SMARCA4 rep3", #/ column title 
                  column_title_gp = gpar(fontsize = 12, fontfamily = "sans"),
                  use_raster = TRUE, raster_quality = 10, raster_device = "png",
                  #/ turn off background colors
                  rect_gp = gpar(col = "transparent"), axis_name_gp = gpar(fontsize = 15),
                  #/ legend options:
                  heatmap_legend_param = list(
                    labels_gp = gpar(font = 15),
                    title = "normalized counts"),
                  #/ options for the profile plot on top of the heatmap:
                  top_annotation = HeatmapAnnotation(enriched = anno_enriched(
                    gp = gpar(col = 1:4, lty = 1:4, lwd=2), axis_param = list(gp = gpar(fontsize =15)),
                    col="black", ylim = c(0,20)))) +
  EnrichedHeatmap(mat = normMatrix_Ig_3, name = "IgG",
                  pos_line = FALSE, #/ no dashed lines around the start
                  border = FALSE,   #/ no box around heatmap
                  col = col_fun_Ig,    #/ color gradients from above
                  column_title = "IgG rep3", #/ column title 
                  column_title_gp = gpar(fontsize = 12, fontfamily = "sans"),
                  use_raster = TRUE, raster_quality = 10, raster_device = "png",
                  #/ turn off background colors
                  rect_gp = gpar(col = "transparent"), axis_name_gp = gpar(fontsize = 15),
                  #/ legend options:
                  heatmap_legend_param = list(
                    labels_gp = gpar(font = 15),
                    title = "normalized counts"),
                  #/ options for the profile plot on top of the heatmap:
                  top_annotation = HeatmapAnnotation(enriched = anno_enriched(
                    gp = gpar(col = 1:4, lty = 1:4, lwd=2), axis_param = list(gp = gpar(fontsize =15)),
                    col="black", ylim = c(0,20))))


pdf("unique_reads/Enriched_heatmaps/PBRM1_non-specific_cen_peaks_rep3_SWI_SNF.pdf", width = 12, height = 12)


draw(EH_NS,                                 #/ plot the heatmap from above 
     heatmap_legend_side = "right",     #/ we want the legend below the heatmap
     annotation_legend_side = "right",  #/ legend on the bottom side
     padding = unit(c(4, 4, 4, 4), "mm"), #/ some padding to avoid labels beyond plot borders
     ht_gap = unit(c(10, 10, 10), "mm"))



dev.off()

tiff('unique_reads/Enriched_heatmaps/PBRM1_non-specific_cen_peaks_rep3_SWI_SNF.tiff', units="in", width=12, height=12, res=300, compression = 'lzw')

draw(EH_NS,                                 #/ plot the heatmap from above 
     heatmap_legend_side = "right",     #/ we want the legend below the heatmap
     annotation_legend_side = "right",  #/ legend on the bottom side
     padding = unit(c(4, 4, 4, 4), "mm"), #/ some padding to avoid labels beyond plot borders
     ht_gap = unit(c(10, 10, 10), "mm"))



dev.off()

emf('unique_reads/Enriched_heatmaps/PBRM1_non-specific_cen_peaks_rep3_SWI_SNF.emf', units="in", width=12, height=12, coordDPI=300)

draw(EH_NS,                                 #/ plot the heatmap from above 
     heatmap_legend_side = "right",     #/ we want the legend below the heatmap
     annotation_legend_side = "right",  #/ legend on the bottom side
     padding = unit(c(4, 4, 4, 4), "mm"), #/ some padding to avoid labels beyond plot borders
     ht_gap = unit(c(10, 10, 10), "mm"))



dev.off()





# KO specific peaks

# make targets file from KO specific peaks - centre around peak - 5kb either side

targets_KO_spec <- makeGRangesFromDataFrame(
  df = fread("unique_reads/MACS2_peaks/overlap_2/CEN_peaks/KO_specific_cen_peaks_merge.bed", header = FALSE, data.table = FALSE),
  seqnames.field = "V1", start.field = "V2", end.field = "V3")

ExtendSize <- 5000
targets_KO_spec.extended  <- resize(targets_KO_spec, fix = "center", width = ExtendSize*2)


# load bigwig files

# CUT&RUN

# WT SMARCA4

BigWig_Sm_3 <- rtracklayer::import("unique_reads/bigwig/WT_SMARCA4_LH_r3.bs10.bw",
                                   format = "BigWig",
                                   selection = BigWigSelection(targets_KO_spec.extended))

normMatrix_Sm_3 <- normalizeToMatrix(signal = BigWig_Sm_3,
                                     target = resize(targets_KO_spec, fix = "center", width = 1),
                                     background = 0,
                                     keep = c(0, 0.99),      #/ minimal value to the 99th percentile
                                     target_ratio = 0,
                                     mean_mode = "coverage",       #/ see ?EnrichedHeatmap on other options
                                     value_column = "score", #/ = the name of the 4th column of the bigwig
                                     extend = ExtendSize)

col_fun_Sm = circlize::colorRamp2(quantile(normMatrix_Sm_3, c(0, .99)), c("white", "#00AD8A"))

# WT PBRM1

BigWig_Pb_3 <- rtracklayer::import("unique_reads/bigwig/WT_PBRM1_LH_r3.bs10.bw",
                                   format = "BigWig",
                                   selection = BigWigSelection(targets_KO_spec.extended))

normMatrix_Pb_3 <- normalizeToMatrix(signal = BigWig_Pb_3,
                                     target = resize(targets_KO_spec, fix = "center", width = 1),
                                     background = 0,
                                     keep = c(0, 0.99),      #/ minimal value to the 99th percentile
                                     target_ratio = 0,
                                     mean_mode = "coverage",       #/ see ?EnrichedHeatmap on other options
                                     value_column = "score", #/ = the name of the 4th column of the bigwig
                                     extend = ExtendSize)

col_fun_Pb = circlize::colorRamp2(quantile(normMatrix_Pb_3, c(0, .99)), c("white", "#009ADE"))

#PBRM1-KO SMARCA4

BigWig_PKO_3 <- rtracklayer::import("unique_reads/bigwig/PBRM1KO_SMARCA4_LH_r3.bs10.bw",
                                    format = "BigWig",
                                    selection = BigWigSelection(targets_KO_spec.extended))

normMatrix_PKO_3 <- normalizeToMatrix(signal = BigWig_PKO_3,
                                      target = resize(targets_KO_spec, fix = "center", width = 1),
                                      background = 0,
                                      keep = c(0, 0.99),      #/ minimal value to the 99th percentile
                                      target_ratio = 0,
                                      mean_mode = "coverage",       #/ see ?EnrichedHeatmap on other options
                                      value_column = "score", #/ = the name of the 4th column of the bigwig
                                      extend = ExtendSize)

col_fun_PKO = circlize::colorRamp2(quantile(normMatrix_PKO_3, c(0, .99)), c("white", "#50A315"))

# WT IgG

BigWig_Ig_3 <- rtracklayer::import("unique_reads/bigwig/WT_IgG_LrmDupH_r3.bs10.bw",
                                   format = "BigWig",
                                   selection = BigWigSelection(targets_KO_spec.extended))

normMatrix_Ig_3 <- normalizeToMatrix(signal = BigWig_Ig_3,
                                     target = resize(targets_KO_spec, fix = "center", width = 1),
                                     background = 0,
                                     keep = c(0, 0.99),      #/ minimal value to the 99th percentile
                                     target_ratio = 0,
                                     mean_mode = "coverage",       #/ see ?EnrichedHeatmap on other options
                                     value_column = "score", #/ = the name of the 4th column of the bigwig
                                     extend = ExtendSize)

col_fun_Ig = circlize::colorRamp2(quantile(normMatrix_Ig_3, c(0, .99)), c("white", "#919191"))


# Prepare heatmap

#### KO specific cen peaks - SWI/SNF signal

set.seed(123)

EH_KS <- 
  EnrichedHeatmap(mat = normMatrix_PKO_3, name = "PBRM1_KO SMARCA4",
                  pos_line = FALSE, #/ no dashed lines around the start
                  border = FALSE,   #/ no box around heatmap
                  col = col_fun_PKO,    #/ color gradients from above
                  column_title = "PBRM1-KO SMARCA4 rep3", #/ column title 
                  column_title_gp = gpar(fontsize = 12, fontfamily = "sans"),
                  use_raster = TRUE, raster_quality = 10, raster_device = "png",
                  #/ turn off background colors
                  rect_gp = gpar(col = "transparent"), axis_name_gp = gpar(fontsize = 15),
                  #/ legend options:
                  heatmap_legend_param = list(
                    labels_gp = gpar(font = 15),
                    title = "normalized counts"),
                  #/ options for the profile plot on top of the heatmap:
                  top_annotation = HeatmapAnnotation(enriched = anno_enriched(
                    gp = gpar(col = 1:4, lty = 1:4, lwd=2), axis_param = list(gp = gpar(fontsize =15)),
                    col="black", ylim = c(0,15)))) +
  EnrichedHeatmap(mat = normMatrix_Sm_3, name = "SMARCA4", 
                  pos_line = FALSE, #/ no dashed lines around the start
                  border = FALSE,   #/ no box around heatmap
                  col = col_fun_Sm,    #/ color gradients from above
                  column_title = "SMARCA4 rep3", #/ column title 
                  column_title_gp = gpar(fontsize = 12, fontfamily = "sans"),
                  use_raster = TRUE, raster_quality = 10, raster_device = "png",
                  #/ turn off background colors
                  rect_gp = gpar(col = "transparent"), axis_name_gp = gpar(fontsize = 15),
                  #/ legend options:
                  heatmap_legend_param = list(
                    labels_gp = gpar(font = 15),
                    title = "normalized counts"),
                  #/ options for the profile plot on top of the heatmap:
                  top_annotation = HeatmapAnnotation(enriched = anno_enriched(
                    gp = gpar(col = 1:4, lty = 1:4, lwd=2), axis_param = list(gp = gpar(fontsize =15)),
                    col="black", ylim = c(0,15)))) +
  EnrichedHeatmap(mat = normMatrix_Pb_3, name = "PBRM1",
                  pos_line = FALSE, #/ no dashed lines around the start
                  border = FALSE,   #/ no box around heatmap
                  col = col_fun_Pb,    #/ color gradients from above
                  column_title = "PBRM1 rep3", #/ column title 
                  column_title_gp = gpar(fontsize = 12, fontfamily = "sans"),
                  use_raster = TRUE, raster_quality = 10, raster_device = "png",
                  #/ turn off background colors
                  rect_gp = gpar(col = "transparent"), axis_name_gp = gpar(fontsize = 15),
                  #/ legend options:
                  heatmap_legend_param = list(
                    labels_gp = gpar(font = 15),
                    title = "normalized counts"),
                  #/ options for the profile plot on top of the heatmap:
                  top_annotation = HeatmapAnnotation(enriched = anno_enriched(
                    gp = gpar(col = 1:4, lty = 1:4, lwd=2), axis_param = list(gp = gpar(fontsize =15)),
                    col="black", ylim = c(0,15)))) +
  EnrichedHeatmap(mat = normMatrix_Ig_3, name = "IgG", 
                  pos_line = FALSE, #/ no dashed lines around the start
                  border = FALSE,   #/ no box around heatmap
                  col = col_fun_Ig,    #/ color gradients from above
                  column_title = "IgG rep3", #/ column title 
                  column_title_gp = gpar(fontsize = 12, fontfamily = "sans"),
                  use_raster = TRUE, raster_quality = 10, raster_device = "png",
                  #/ turn off background colors
                  rect_gp = gpar(col = "transparent"), axis_name_gp = gpar(fontsize = 15),
                  #/ legend options:
                  heatmap_legend_param = list(
                    labels_gp = gpar(font = 15),
                    title = "normalized counts"),
                  #/ options for the profile plot on top of the heatmap:
                  top_annotation = HeatmapAnnotation(enriched = anno_enriched(
                    gp = gpar(col = 1:4, lty = 1:4, lwd=2), axis_param = list(gp = gpar(fontsize =15)),
                    col="black", ylim = c(0,15))))


pdf("unique_reads/Enriched_heatmaps/KO_specific_cen_peaks_rep3_SWI_SNF.pdf", width = 12, height = 12)


draw(EH_KS,                                 #/ plot the heatmap from above 
     heatmap_legend_side = "right",     #/ we want the legend below the heatmap
     annotation_legend_side = "right",  #/ legend on the bottom side
     padding = unit(c(4, 4, 4, 4), "mm"), #/ some padding to avoid labels beyond plot borders
     ht_gap = unit(c(10, 10, 10), "mm"))



dev.off()

tiff('unique_reads/Enriched_heatmaps/KO_specific_cen_peaks_rep3_SWI_SNF.tiff', units="in", width=12, height=12, res=300, compression = 'lzw')

draw(EH_KS,                                 #/ plot the heatmap from above 
     heatmap_legend_side = "right",     #/ we want the legend below the heatmap
     annotation_legend_side = "right",  #/ legend on the bottom side
     padding = unit(c(4, 4, 4, 4), "mm"), #/ some padding to avoid labels beyond plot borders
     ht_gap = unit(c(10, 10, 10), "mm"))



dev.off()

emf('unique_reads/Enriched_heatmaps/KO_specific_cen_peaks_rep3_SWI_SNF.emf', units="in", width=12, height=12, coordDPI=300)

draw(EH_KS,                                 #/ plot the heatmap from above 
     heatmap_legend_side = "right",     #/ we want the legend below the heatmap
     annotation_legend_side = "right",  #/ legend on the bottom side
     padding = unit(c(4, 4, 4, 4), "mm"), #/ some padding to avoid labels beyond plot borders
     ht_gap = unit(c(10, 10, 10), "mm"))



dev.off()

