library(QDNAseq)
library(dplyr)
options(future.globals.maxSize = 1048576000)
future::plan("multiprocess")

# ChIp-seq Metadata
metadata_chipseq = data.frame(
  CellType = c("GICAN", "SH-EP", "SK-N-AS", "GIMEN", "SK-N-SH", "NB69", "SJNB12",
               "SH-SY5Y", "SJNB1", "SK-N-FI", "CLB-GA", "NB-EBc1",
               "LAN1", "CLB-PE", "SK-N-DZ", "CLB-CAR", "CLB-MA",
               "IMR32", "CHP212", "SJNB8", "TR14", "SK-N-BE2-C",
               "N206", "SJNB6", "CLB-BER-Lud")
  
)
metadata_chipseq$Class = NA
metadata_chipseq = 
  metadata_chipseq %>%
  mutate(Class = 
           ifelse(CellType %in% c("GICAN", "SH-EP", "SK-N-AS", "GIMEN", "SK-N-SH", "NB69", "SJNB12", "SJNB12"), 
                  "noMYCN", 
                  Class)) %>%
  mutate(Class = 
           ifelse(CellType %in% c("SH-SY5Y", "SJNB1", "SK-N-FI", "CLB-GA", "NB-EBc1"), 
                  "lowMYCN", 
                  Class)) %>%
  mutate(Class = 
           ifelse(CellType %in% c("LAN1", "CLB-PE", "SK-N-DZ", "CLB-CAR", "CLB-MA",
                                  "IMR32", "CHP212", "SJNB8", "TR14", "SK-N-BE2-C",
                                  "N206", "SJNB6", "CLB-BER-Lud"), 
                  "MNA", 
                  Class))

metadata_chipseq$bam_fname = paste0("/Volumes/Elements/nb-cl-chipseq-results/bam/Boeva_", metadata_chipseq$CellType, "_Input.trimmed.bwa_hg19.rmdup.bam")
metadata_chipseq$output_bed_fname = paste0("/Volumes/Elements/nb-cl-chipseq-qdnaseq/Boeva_", metadata_chipseq$CellType, "_Input.trimmed.bwa_hg19.rmdup.bam.qdnaseq.bed")
metadata_chipseq$output_pdf_fname = paste0("/Volumes/Elements/nb-cl-chipseq-qdnaseq/Boeva_", metadata_chipseq$CellType, "_Input.trimmed.bwa_hg19.rmdup.bam.qdnaseq.pdf")

# this is just a very special case
metadata_chipseq[11, "bam_fname"] = paste0("/Volumes/Elements/nb-cl-chipseq-results/bam/Boeva_CLB-GA_rep2_Input.trimmed.bwa_hg19.rmdup.bam")

# make sure all files exist
sum(!file.exists(metadata_chipseq$bam_fname))

# Finished i=10. Afterewards broke at CLB-GA, bc did not find the bam file.
# 17:16 Finished i=14. 
# 18:39 Finished i=16
# 20:05 Finished i=18
# 23:13 Finished i=22


for (i in 15:nrow(metadata_chipseq)){
  bins <- getBinAnnotations(binSize = 1, genome="hg19")
  r <- binReadCounts(bins, bamfiles = metadata_chipseq[i, "bam_fname"], 
                     minMapq=20, isSecondaryAlignment=FALSE)
  r <- applyFilters(r, residual = F, blacklist = T, mappability = F, 
                    bases = F, chromosomes = c("X", "Y", "MT"))
  r <- estimateCorrection(r)
  r <- correctBins(r)
  r <- normalizeBins(r)
  r <- segmentBins(r, alpha = 0.01, transformFun = "sqrt")
  
  # Save results to bed file
  exportBins(r, file=metadata_chipseq[i, "output_bed_fname"], 
             format="bed", type = "segments", 
             filter = T, logTransform = F, digits = 2)
  
  # Plot bins and segmentation
  # if (file.exists(as.character(metadata_chipseq[i, "output_pdf_fname"])) )file.remove(as.character(metadata_chipseq[i, "output_pdf_fname"]))
  # pdf(file = as.character(metadata_chipseq[i, "output_pdf_fname"]))
  # plot(r)
  # dev.off()
  
  print("Finished:")
  print(i)
}
