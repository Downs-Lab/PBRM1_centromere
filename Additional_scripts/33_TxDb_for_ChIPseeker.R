setwd(system("pwd", intern = T) )


if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos = "http://cran.us.r-project.org")

BiocManager::install(c("txdbmaker", "rtracklayer", "GenomeInfoDb", "AnnotationDbi"))



library(rtracklayer)
library(GenomeInfoDb)
library(txdbmaker)
library(AnnotationDbi)

gff <- import("genomes/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz")

seqlevels(gff)


chrominfo <- getChromInfoFromNCBI("T2T-CHM13v2.0")
seqlevels(gff) <- setNames(chrominfo$SequenceName, chrominfo$RefSeqAccn)
seqlevels(gff)

seqinfo(gff) <- Seqinfo(genome="T2T-CHM13v2.0")
seqinfo(gff)

seqlevelsStyle(gff) <- "UCSC"

# Make TxDb file from annotation files for use in ChIPseeker
## This will emit 3 warnings that can be ignored.
txdb <- makeTxDbFromGRanges(gff, taxonomyId=9606)

txdb

saveDb(txdb, "genomes/txdb_t2t")

# now you can reload the database without redoing earlier steps

txdb <- loadDb("genomes/txdb_t2t")

