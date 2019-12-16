library(DESeq2)
library(ggplot2)
library(dplyr)
library(data.table)

metadata = read.table("/Volumes/Elements/MYCNAmplicon/Data/Boeva_RNAseq_Metadata.csv", header=T, sep=";")

l <- lapply(as.character(metadata[, "featureCounts_fname"]), function(fn) read.table(file.path(fn), skip=2))
if (!all(sapply(l, function(a) all(a$V1 == l[[1]]$V1)))) stop("Gene IDs (first column) differ between files.")
tbl <- sapply(l, function(a) a$V7)
colnames(tbl) <- metadata[, "Sample"]
rownames(tbl) <- l[[1]]$V1
rownames(metadata) <- metadata$Sample
dds <- 
  DESeqDataSetFromMatrix(
    countData = tbl,
    colData = metadata[,c("Sample", "CellType", "Study", "MYCNStatus", "Library")],
    design = ~ CellType)
dds <- dds[ rowSums(counts(dds)) > 1,]

dds <- estimateSizeFactors(dds)
counts.norm = counts(dds, normalized=TRUE)

counts.norm %>% 
  write.table(paste0("/Volumes/Elements/MYCNAmplicon/Workspace/boeva_rnaseq.sizeFactorNormalized.txt"), 
              sep="\t", col.names=T, row.names=T, quote=F)

vsd <- vst(dds, blind=FALSE)

save.image(paste0("/Volumes/Elements/MYCNAmplicon/Workspace/boeva_rnaseq.Rdata"))
