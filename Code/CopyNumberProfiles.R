rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggforce)
source("/Volumes/Transcend/PalmTrees/Analysis/Code/hasOverlap.R")
source("/Volumes/Transcend/PalmTrees/Analysis/Code/CustomThemes.R")
source("/Volumes/Transcend/PalmTrees/Analysis/Code/ggGrandLinear.R")
load("/Volumes/Transcend/PalmTrees/Analysis/WorkspaceData/ascat_cnv.Rdata")

clinical_data = read.table("/Volumes/Elements/MYCNAmplicon/Data/ClinicalData.csv",
                           header=T, sep=",")
mna_samples = clinical_data[clinical_data$Risk == "MNA","Sample"] %>% as.character() %>% unique()

ascat_cnv =
  ascat_cnv %>% filter(Sample %in% mna_samples)

samples = unique(ascat_cnv$Sample)
for (i in 1:length(samples)){
  this_sample = samples[i]
  this_ascat_cnv_gr = ascat_cnv_gr[ascat_cnv_gr$Sample == this_sample]
  this_ascat_cnv = data.frame(this_ascat_cnv_gr) %>% as_tibble()
  
  this_ascat_cnv %>%
    filter(seqnames == "chr2") %>%
    ggplot(aes(x = start, y = TumorTotalCopyNumber)) +
    geom_segment(aes(xend = end, yend = TumorTotalCopyNumber), size=3) +
    xlab("chr2") + ylab("Copy Number (ASCAT)") +
    scale_y_log10() +
    annotate(xmin=16080683, xmax=16087129, ymin=0, ymax=Inf, geom="rect", color=NA, fill="black") + # MYCN
    annotate(xmin=16375276, xmax=16385344, ymin=0, ymax=Inf, geom="rect", color=NA, fill="red") + # MYCN
    facet_zoom(xlim=c(10000000,20000000), zoom.size=1) + 
    theme_kons1() +
    theme(strip.background = element_rect(fill="grey90", color=NA)) +
    ggsave(paste0("/Volumes/Elements/MYCNAmplicon/Figures/ASCATFigures/", this_sample, "_ASCAT_chr2.pdf"),
           height = 5, width = 10, useDingbats = F)
    
  ggGrandLinearBed(this_ascat_cnv$seqnames, this_ascat_cnv$start, this_ascat_cnv$end, this_ascat_cnv$TumorTotalCopyNumber) +
      theme_kons1() +
      ylab("Copy Number (ASCAT)") +
      ggtitle(this_sample) +
    ggsave(paste0("/Volumes/Elements/MYCNAmplicon/Figures/ASCATFigures/", this_sample, "_ASCAT_GrandLinear.pdf"),
           height = 3, width = 10, useDingbats = F)
  
    # freec %>%
    # as_tibble() %>%
    # mutate(type = factor(type, levels = c("loss", "normal", "gain"))) %>% 
    # filter(seqnames == "chr2", ) %>%
    # filter(copy.number < 600) %>% 
    #  +
    # ggtitle(sample_name) +
    # 
    # theme_kons1() + 
    # scale_color_manual(values=unlist(list("loss" = "steelblue", "normal" = "grey75", "gain" = "firebrick2"))) + 
    # guides(color=F) +
    # ggsave(paste0(freec_fname, "chr2.pdf"), height=5, width=5)  
  
}


