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

gene_of_interest = "CDK4"
gene_of_interest_ensembl = "ENSG00000135446"
gene_of_interest_chr = "chr12"
gene_of_interest_start = 58141510
gene_of_interest_end = 58149796
window_of_interest = 500000

cdk4_amplified_profiles = 
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
seqlevels(cdk4_amplified_profiles, pruning.mode = "coarse") = seqlevels(BSgenome.Hsapiens.UCSC.hg19)

# CDK4 chr12:58,141,510-58,149,796
cdk4_granges = GRanges(
  seqnames = c("chr12"),
  ranges = IRanges(
    start = c(58141510+2000),
    end = c(58149796-2000)
  )
)
arrayDataIsCDK4ampl = vector()
for (sample_idx in 1:length(cdk4_amplified_profiles)){
  if (length(findOverlaps(cdk4_granges, cdk4_amplified_profiles[sample_idx])) > 0) arrayDataIsCDK4ampl = c(arrayDataIsCDK4ampl, sample_idx)
}
cdk4_amplified_profiles = cdk4_amplified_profiles[arrayDataIsCDK4ampl]

# how many samples do we have?
length(cdk4_amplified_profiles)

load("/Volumes/Transcend/PalmTrees/Analysis/WorkspaceData/CNVData.Rdata")
wgs_cdk4_amplicon_samples = 
  cnv_gr %>% 
  filter(CopyNumber > 9) %>%
  filter_by_overlaps(cdk4_granges) %>%
  .$Sample %>% unique
wgs_cdk4_profiles = 
  cnv_gr %>%
  filter(Sample %in% wgs_cdk4_amplicon_samples)

cdk4_amplified_profiles_df = 
  cdk4_amplified_profiles %>% 
  unlist() %>% as_tibble()
cdk4_amplified_profiles_df$Sample = rownames(cdk4_amplified_profiles)

wgs_cdk4_profiles = GRangesList()
for (sample in wgs_cdk4_amplicon_samples){
  wgs_cdk4_profiles[[sample]] = GRanges(cnv_gr %>% filter(Sample == sample, CopyNumber > 9))
}

cdk4_amplified_profiles = c(cdk4_amplified_profiles, wgs_cdk4_profiles)

bins = tileGenome(seqlengths = seqlengths(BSgenome.Hsapiens.UCSC.hg19), 
                  tilewidth = 10000,
                  cut.last.tile.in.chrom = T)
bins_focus = bins %>% plyranges::filter(start > 55000000, end < 65000000)
cdk4_amplified_profiles_cov = lapply(cdk4_amplified_profiles, coverage)
cdk4_amplified_profiles_binned = lapply(cdk4_amplified_profiles_cov, function (this_profile_cov) binnedAverage(bins, numvar = this_profile_cov, varname = "isAmp"))
this_profile = cdk4_amplified_profiles_binned[[1]]
for (i in 1:length(cdk4_amplified_profiles_binned)){
  cdk4_amplified_profiles_binned[[i]]$Name = names(cdk4_amplified_profiles_binned)[i]
}
cdk4_amplified_profiles_binned_tb = do.call(rbind, lapply(cdk4_amplified_profiles_binned, as_tibble))

samples = samples %>% filter(Name %in% unique(cdk4_amplified_profiles_binned_tb$Name))
mna_profiles_binned_tb = 
  mna_profiles_binned_tb %>%
  full_join(samples, by= "Name")

cdk4_amplified_profiles_binned_tb %>%
  filter(seqnames == "chr12") %>% 
  group_by(start) %>%
  summarise(PercentSamples = 100*mean(isAmp),
            nSamples = sum(isAmp)) %>%
  ungroup() %>%
  ggplot(aes(x=start, y=PercentSamples)) +
  geom_line() + 
  facet_zoom(xlim = c(55000000, 65000000)) +
  theme_kons1() + 
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ylab("Samples [%]") + 
  xlab("chr2") + 
  annotate(geom="rect", xmin=gene_of_interest_start, xmax=gene_of_interest_end, ymin=0, ymax=Inf, color=NA, fill="firebrick3", alpha=1) +
  ggtitle("Amplified Regions around CDK4 (Depuydt et al. 2018 + Berlin/Peifer Cohort)") + 
  ggsave("/Volumes/Elements/Celine_CDK4/CDK4Amplicon_Consensus.pdf",
         height = 5, width = 5, useDingbats=F)


# ------------------------------------------------------------------------------
# MDM2
# ------------------------------------------------------------------------------

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
arrayDataIsMDM2ampl = vector()
for (sample_idx in 1:length(mdm2_amplified_profiles)){
  if (length(findOverlaps(mdm2_granges, mdm2_amplified_profiles[sample_idx])) > 0) arrayDataIsMDM2ampl = c(arrayDataIsMDM2ampl, sample_idx)
}
mdm2_amplified_profiles = mdm2_amplified_profiles[arrayDataIsMDM2ampl]

# how many samples do we have?
length(mdm2_amplified_profiles)


