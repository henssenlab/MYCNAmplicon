rm(list=ls())

library(plyranges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(biomaRt)
library(ggforce)
source("/Volumes/Elements/MYCNAmplicon/Code/CustomThemes.R")

setwd("~/Desktop/copynumber_HR_NB-master/")
samples <- read.delim("samples_dd.txt")
samples$Name = paste0("Sample", as.character(samples$Name))
colnames(samples)[colnames(samples) == "MYCN"] = "MYCNStatus"
profiles_df <- read.delim("profiles_dd.txt")
profiles_df$Name = paste0("Sample", as.character(profiles_df$Name))

gene_of_interest = "MDM2"
gene_of_interest_ensembl = "ENSG00000135679"
gene_of_interest_chr = "chr12"
gene_of_interest_start = 69201956
gene_of_interest_end = 69239214
window_of_interest = 500000

mdm2_amplified_profiles = 
  profiles_df %>%
  full_join(samples, by="Name") %>%
  mutate(chromosome = paste0("chr", as.character(chromosome))) %>%
  mutate(chromosome = as.character(chromosome)) %>% 
  filter(annotation == "ampl") %>%
  makeGRangesListFromDataFrame(keep.extra.columns = T,
                               split.field = "Name",
                               seqnames.field = "chromosome",
                               start.field = "min",
                               end.field = "max")
seqlevels(mdm2_amplified_profiles, pruning.mode = "coarse") = seqlevels(BSgenome.Hsapiens.UCSC.hg19)

mdm2_granges = GRanges(
  seqnames = c("chr12"),
  ranges = IRanges(
    start = c(69201956),
    end = c(69239214)
  )
)
arrayDataIsmdm2ampl = vector()
for (sample_idx in 1:length(mdm2_amplified_profiles)){
  if (length(findOverlaps(mdm2_granges, mdm2_amplified_profiles[sample_idx])) > 0) arrayDataIsmdm2ampl = c(arrayDataIsmdm2ampl, sample_idx)
}
mdm2_amplified_profiles = mdm2_amplified_profiles[arrayDataIsmdm2ampl]

# how many samples do we have?
length(mdm2_amplified_profiles)

load("/Volumes/Transcend/PalmTrees/Analysis/WorkspaceData/CNVData.Rdata")
wgs_mdm2_amplicon_samples = 
  cnv_gr %>% 
  filter(CopyNumber > 9) %>%
  filter_by_overlaps(mdm2_granges) %>%
  .$Sample %>% unique
wgs_mdm2_profiles = 
  cnv_gr %>%
  filter(Sample %in% wgs_mdm2_amplicon_samples)

mdm2_amplified_profiles_df = 
  mdm2_amplified_profiles %>% 
  unlist() %>% as_tibble()
mdm2_amplified_profiles_df$Sample = rownames(mdm2_amplified_profiles)

wgs_mdm2_profiles = GRangesList()
for (sample in wgs_mdm2_amplicon_samples){
  wgs_mdm2_profiles[[sample]] = GRanges(cnv_gr %>% filter(Sample == sample, CopyNumber > 9))
}

mdm2_amplified_profiles = c(mdm2_amplified_profiles, wgs_mdm2_profiles)

bins = tileGenome(seqlengths = seqlengths(BSgenome.Hsapiens.UCSC.hg19), 
                  tilewidth = 10000,
                  cut.last.tile.in.chrom = T)
bins_focus = bins %>% plyranges::filter(start > 65000000, end < 75000000)
mdm2_amplified_profiles_cov = lapply(mdm2_amplified_profiles, coverage)
mdm2_amplified_profiles_binned = lapply(mdm2_amplified_profiles_cov, function (this_profile_cov) binnedAverage(bins, numvar = this_profile_cov, varname = "isAmp"))
this_profile = mdm2_amplified_profiles_binned[[1]]
for (i in 1:length(mdm2_amplified_profiles_binned)){
  mdm2_amplified_profiles_binned[[i]]$Name = names(mdm2_amplified_profiles_binned)[i]
}
mdm2_amplified_profiles_binned_tb = do.call(rbind, lapply(mdm2_amplified_profiles_binned, as_tibble))

samples = samples %>% filter(Name %in% unique(mdm2_amplified_profiles_binned_tb$Name))

mdm2_amplified_profiles_binned_tb %>%
  filter(seqnames == "chr12") %>% 
  group_by(start) %>%
  summarise(PercentSamples = 100*mean(isAmp),
            nSamples = sum(isAmp)) %>%
  ungroup() %>%
  ggplot(aes(x=start, y=PercentSamples)) +
  geom_line() + 
  facet_zoom(xlim = c(65000000, 75000000)) +
  theme_kons1() + 
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ylab("Samples [%]") + 
  xlab("chr2") + 
  annotate(geom="rect", xmin=gene_of_interest_start, xmax=gene_of_interest_end, ymin=0, ymax=Inf, color=NA, fill="firebrick3", alpha=1) +
  ggtitle("Amplified Regions around MDM2 (Depuydt et al. 2018 + Berlin/Peifer Cohort)") + 
  ggsave("/Volumes/Elements/Celine_CDK4/MDM2_Amplicon_Consensus.pdf",
         height = 5, width = 5, useDingbats=F)


