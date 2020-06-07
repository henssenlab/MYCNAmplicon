library(plyranges)
library(ggplot2)
library(dplyr)
source("/Volumes/Elements/MYCNAmplicon/Code/CustomThemes.R")

metadata = data.frame(
  CellLine = c("CHP-212", "IMR-5-75", "Kelly", "NGP"),
  BigWig = c("/Volumes/Elements/CHP-212.ONT-flye-assembly-meta-genomesize1g/reads2assembly.ngmlr.cov.bw",
             "/Volumes/Elements/IMR-5-75.ONT-flye-assembly-meta-genomesize1g/reads2assembly.ngmlr.cov.bw",
             "/Volumes/Elements/Kelly.ONT-flye-assembly-meta-genomesize1g/reads2assembly.ngmlr.cov.bw",
             "/Volumes/Elements/NGP.ONT-flye-assembly-meta-genomesize1g/reads2assembly.ngmlr.cov.bw"),
  MYCNContig = c("contig_1695",
                 "contig_499",
                 "contig_31",
                 "contig_45")
)
metadata$CellLine = as.character(metadata$CellLine)
metadata$BigWig = as.character(metadata$BigWig)
metadata$MYCNContig = as.character(metadata$MYCNContig)


for (i in 1:nrow(metadata)){
  bw = read_bigwig(metadata[i,"BigWig"])
  bw %>% 
    plyranges::filter(seqnames == metadata[i,"MYCNContig"]) %>%
    as_tibble() %>%
    ggplot(aes(x=start, y=score)) +
    geom_line() +
    ggtitle(paste0(metadata[i, "CellLine"], "\nMYCN-containing contig (", metadata[i, "MYCNContig"], ")")) +
    expand_limits(x = 0, y = 0) +
    xlab("Position in contig [bp]") +
    ylab("Read Coverage") +
    #ylim(0,400) +
    theme_kons2() +
    ggsave(paste0("/Volumes/Elements/MYCNAmplicon/Results/ContigCoverage_",
                  metadata[i, "CellLine"],
                  "_",
                  metadata[i,"MYCNContig"], 
                  ".pdf"),
           height = 2, width = 3) 
}
