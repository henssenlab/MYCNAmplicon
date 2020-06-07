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

gene_of_interest = "MYCN"
gene_of_interest_ensembl = "ENSG00000134323"
gene_of_interest_chr = "chr2"
gene_of_interest_start = 16080683
gene_of_interest_end = 16087129
window_of_interest = 500000

# ------------------------------------------------------------------------------
# Take only amplified segments on chr2 for clinically annotated MNA samples 
# ------------------------------------------------------------------------------

mna_profiles = profiles_df %>%
  full_join(samples, by="Name") %>%
  filter(MYCNStatus == 1) %>%
  mutate(chromosome = paste0("chr", as.character(chromosome))) %>%
  filter(chromosome == "chr2") %>%
  mutate(chromosome = as.character(chromosome)) %>% 
  filter(annotation == "ampl") %>%
  makeGRangesListFromDataFrame(keep.extra.columns = T,
                               split.field = "Name",
                               seqnames.field = "chromosome",
                               start.field = "min",
                               end.field = "max", 
                               seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg19))
seqlevels(mna_profiles, pruning.mode = "coarse") = seqlevels(BSgenome.Hsapiens.UCSC.hg19)[2]


# ------------------------------------------------------------------------------
# Filter out samples that are clincially annotated as MNA but do not show MYCN
# amplification in the array data
# ------------------------------------------------------------------------------

mycn_granges = GRanges(
  seqnames = c("chr2"),
  ranges = IRanges(
    start = c(16080683),
    end = c(16087129)
  )
)
arrayDataIsMNA = vector()
for (sample_idx in 1:length(mna_profiles)){
  if (length(findOverlaps(mycn_granges, mna_profiles[sample_idx])) > 0) arrayDataIsMNA = c(arrayDataIsMNA, sample_idx)
}
mna_profiles = mna_profiles[arrayDataIsMNA]

# how many samples do we have?
length(mna_profiles)

mna_profiles %>%
  write_bed("/Volumes/Elements/MYCNAmplicon/Results/MNAProfilesChr2.bed")

# ------------------------------------------------------------------------------
# Bin amplification data for MNA samples
# ------------------------------------------------------------------------------

bins = tileGenome(seqlengths = seqlengths(BSgenome.Hsapiens.UCSC.hg19)[2], 
                  tilewidth = 10000,
                  cut.last.tile.in.chrom = T)
bins_focus = bins %>% plyranges::filter(start > 5000000, end < 35000000)
mna_profiles_cov = lapply(mna_profiles, coverage)
mna_profiles_binned = lapply(mna_profiles_cov, function (this_profile_cov) binnedAverage(bins, numvar = this_profile_cov, varname = "isAmp"))
this_profile = mna_profiles_binned[[1]]

for (i in 1:length(mna_profiles_binned)){
  mna_profiles_binned[[i]]$Name = names(mna_profiles_binned)[i]
}
mna_profiles_binned_tb = do.call(rbind, lapply(mna_profiles_binned, as_tibble))


# ------------------------------------------------------------------------------
# Plot amplification meta-profile
# ------------------------------------------------------------------------------

# How large is the region above 90%?
mna_profiles_binned_tb %>%
  group_by(start) %>%
  summarise(PercentSamples = 100*mean(isAmp),
            nSamples = sum(isAmp)) %>%
  ungroup() %>% 
  filter(PercentSamples > 90) %>%
  View
# last bin with >90% --> 16,370,001-16,380,000
# MYCN = 16080683-16087129
# distance from MYCN gene end: 16380000 - 16087129 = 292,871

mna_profiles_binned_tb %>%
  group_by(start) %>%
  summarise(PercentSamples = 100*mean(isAmp),
            nSamples = sum(isAmp)) %>%
  ungroup() %>%
  ggplot(aes(x=start, y=PercentSamples)) +
  geom_line() + 
  facet_zoom(xlim = c(15000000, 17000000)) +
  theme_kons1() + 
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ylab("Samples [%]") + 
  xlab("chr2") + 
  ggtitle("Amplified Regions on chr2 (Depuydt et al. 2018)") +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/AmplifiedOnChr2_Zoom_All.pdf",
         height = 5, width = 5, useDingbats=F)

samples = samples %>% filter(Name %in% unique(mna_profiles_binned_tb$Name))
mna_profiles_binned_tb = mna_profiles_binned_tb %>%
  full_join(samples, by= "Name") 

# copy number profile by age quartile
mna_profiles_binned_tb %>%
  filter(!is.na(Age)) %>% 
  mutate(AgeQuartile = as.factor(ntile(Age, 4))) %>%
  group_by(start, AgeQuartile) %>%
  summarise(PercentSamples = 100*mean(isAmp),
            nSamples = sum(isAmp)) %>%
  ungroup() %>%
  ggplot(aes(x=start, y=PercentSamples, color=AgeQuartile)) +
  geom_line() + 
  xlim(gene_of_interest_start-500000,gene_of_interest_end + 500000) +
  theme_kons2() + 
  scale_color_brewer(palette = "Set1") +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ylab("Samples [%]") +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/ConsensusAmplicon_byAge.pdf",
         height = 2, width = 2, useDingbats=F)
mna_profiles_binned_tb %>%
  filter(!is.na(Age)) %>% 
  mutate(AgeQuartile = as.factor(ntile(Age, 4))) %>%
  group_by(AgeQuartile) %>%
  summarise(n=n_distinct(Name))
mna_profiles_binned_tb %>%
  filter(!is.na(Age)) %>% 
  mutate(AgeQuartile = as.factor(ntile(Age, 4))) %>%
  group_by(start, AgeQuartile) %>%
  summarise(PercentSamples = 100*mean(isAmp),
            nSamples = sum(isAmp)) %>%
  ungroup() %>%
  ggplot(aes(x=start, y=PercentSamples, color=AgeQuartile)) +
  geom_line() + 
  xlim(gene_of_interest_start-2000000,gene_of_interest_end + 2000000) +
  theme_kons2() + 
  scale_color_brewer(palette = "Set1") +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ylab("Samples [%]") +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/ConsensusAmplicon_byAge_ZoomOut.pdf",
         height = 2, width = 3, useDingbats=F)

# Profile by long-term survival vs. not
mna_profiles_binned_tb %>%
  filter(!is.na(OS)) %>% 
  mutate(longTermSurvival = (OS == 0) & (OStime > 5*365)) %>% 
  group_by(start, longTermSurvival) %>%
  summarise(PercentSamples = 100*mean(isAmp),
            nSamples = sum(isAmp)) %>%
  ungroup() %>%
  ggplot(aes(x=start, y=PercentSamples, color=longTermSurvival)) +
  geom_line() + 
  xlim(gene_of_interest_start-500000,gene_of_interest_end + 500000) +
  theme_kons2() + 
  scale_color_brewer(palette = "Set1") +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ylab("Samples [%]") +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/ConsensusAmplicon_byLongTermSurvival.pdf",
         height = 2, width = 2, useDingbats=F)
mna_profiles_binned_tb %>%
  filter(!is.na(OS)) %>% 
  mutate(longTermSurvival = (OS == 0) & (OStime > 5*365)) %>% 
  dplyr::select(Name, longTermSurvival) %>%
  group_by(longTermSurvival) %>%
  summarise(n=dplyr::n_distinct(Name))

# Profile by copy number method / platform
mna_profiles_binned_tb %>%
  filter(!is.na(Platform)) %>% 
  group_by(start, Platform) %>%
  summarise(PercentSamples = 100*mean(isAmp),
            nSamples = sum(isAmp)) %>%
  ungroup() %>%
  ggplot(aes(x=start, y=PercentSamples, color=Platform)) +
  geom_line() + 
  xlim(gene_of_interest_start-500000,gene_of_interest_end + 500000) +
  theme_kons2() + 
  scale_color_brewer(palette = "Set1") +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ylab("Samples [%]") +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/ConsensusAmplicon_byMethod.pdf",
         height = 2, width = 2, useDingbats=F)
mna_profiles_binned_tb %>%
  filter(!is.na(Platform)) %>% 
  dplyr::select(Name, Platform) %>%
  group_by(Platform) %>%
  summarise(n=dplyr::n_distinct(Name))

# Profile by copy number method / platform
mna_profiles_binned_tb %>%
  mutate(chr17q = as.factor(chr17q)) %>% 
  filter(!is.na(chr17q)) %>% 
  group_by(start, chr17q) %>%
  summarise(PercentSamples = 100*mean(isAmp),
            nSamples = sum(isAmp)) %>%
  ungroup() %>%
  ggplot(aes(x=start, y=PercentSamples, color=chr17q)) +
  geom_line() + 
  xlim(gene_of_interest_start-500000,gene_of_interest_end + 500000) +
  theme_kons2() + 
  scale_color_brewer(palette = "Set1") +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ylab("Samples [%]") +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/ConsensusAmplicon_bychr17q.pdf",
         height = 2, width = 2, useDingbats=F)
mna_profiles_binned_tb %>%
  filter(!is.na(chr17q)) %>% 
  dplyr::select(Name, chr17q) %>%
  group_by(chr17q) %>%
  summarise(n=dplyr::n_distinct(Name))

# Profile by copy number method / platform
mna_profiles_binned_tb %>%
  mutate(chr1p = as.factor(chr1p)) %>% 
  filter(!is.na(chr1p)) %>% 
  group_by(start, chr1p) %>%
  summarise(PercentSamples = 100*mean(isAmp),
            nSamples = sum(isAmp)) %>%
  ungroup() %>%
  ggplot(aes(x=start, y=PercentSamples, color=chr1p)) +
  geom_line() + 
  xlim(gene_of_interest_start-500000,gene_of_interest_end + 500000) +
  theme_kons2() + 
  scale_color_brewer(palette = "Set1") +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ylab("Samples [%]") +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/ConsensusAmplicon_bychr1p.pdf",
         height = 2, width = 2, useDingbats=F)
mna_profiles_binned_tb %>%
  filter(!is.na(chr1p)) %>% 
  dplyr::select(Name, chr1p) %>%
  group_by(chr1p) %>%
  summarise(n=dplyr::n_distinct(Name))

mna_profiles_binned_tb %>%
  group_by(start) %>%
  summarise(PercentSamples = 100*mean(isAmp),
            nSamples = sum(isAmp)) %>%
  ungroup() %>%
  #filter(start>5000000, start<50000000) %>% 
  ggplot(aes(x=start, y=PercentSamples)) +
  geom_line() + 
  xlim(15000000, 17000000) +
  theme_kons1() + 
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ylab("Samples [%]") + 
  xlab("chr2") + 
  ggtitle("Amplified Regions on chr2 (Depuydt et al. 2018)") +
  ggsave("~/Desktop/copynumber_HR_NB-master/AmplifiedOnChr2_MYCNFocus_All.pdf",
         height = 3, width = 5, useDingbats=F)

mna_profiles_binned_tb %>%
  group_by(start) %>%
  summarise(PercentSamples = 100*mean(isAmp),
            nSamples = sum(isAmp)) %>%
  ungroup() %>%
  ggplot(aes(x=start, y=PercentSamples)) +
  geom_line() + 
  xlim(15000000, 17000000) +
  theme_kons1() + 
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ylab("Samples [%]") + 
  xlab("chr2") + 
  ggtitle("Amplified Regions on chr2 (Depuydt et al. 2018)") +
  ggsave("~/Desktop/copynumber_HR_NB-master/AmplifiedOnChr2_MYCNFocus_All.pdf",
         height = 3, width = 5, useDingbats=F)

mna_profiles_binned_tb %>%
  inner_join(samples %>% dplyr::select(Name, Platform), by="Name") %>% 
  group_by(start, Platform) %>%
  summarise(PercentSamples = 100*mean(isAmp),
            nSamples = sum(isAmp)) %>%
  ungroup() %>%
  ggplot(aes(x=start, y=PercentSamples, color=Platform)) +
  scale_color_brewer(palette = "Set1") +
  geom_line() + 
  xlim(15000000, 17000000) +
  theme_kons1() + 
  ylab("Samples [%]") + 
  xlab("chr2") + 
  ggtitle("Amplified Regions on chr2 (Depuydt et al. 2018)") +
  ggsave("~/Desktop/copynumber_HR_NB-master/AmplifiedOnChr2_MYCNFocus_ByPlatform.pdf",
         height = 3, width = 5, useDingbats=F)

mna_profiles_binned_tb %>%
  group_by(start) %>%
  summarise(PercentSamples = 100*mean(isAmp),
            nSamples = sum(isAmp)) %>%
  ungroup() %>%
  ggplot(aes(x=start, y=PercentSamples)) +
  geom_line() + 
  facet_zoom(xlim = c(28000000, 33000000)) +
  theme_kons1() + ylim(0,15) +
  theme(strip.background = element_rect(fill="grey90", color=NA))+
  annotate(xmin=29415640, xmax=30144477, ymin=0, ymax=Inf, geom="rect", color=NA, fill="steelblue", alpha=0.5)+
  ylab("Samples [%]") + 
  xlab("chr2") + 
  ggtitle("Amplified Regions on chr2 (Depuydt et al. 2018)\nFocus on ALK") +
  ggsave("~/Desktop/copynumber_HR_NB-master/AmplifiedOnChr2_ALKFocus_ByPlatform.pdf",
         height = 3, width = 5, useDingbats=F)

# ------------------------------------------------------------------------------
# Plot of copy number meta-profile with gene annotations
# ------------------------------------------------------------------------------

# NBAS chr2:15,307,032-15,701,454
# DDX1 chr2:15,731,302-15,771,235
# MYCN chr2:16,080,683-16,087,129
# GACAT3 chr2:16,190,549-16,225,923
# FAM49A chr2:16,730,727-16,847,599

mna_samples = 
  mna_profiles_binned_tb %>%
  filter(start>16080683, start<16080683+10000, isAmp == 1) %>% 
  .$Name

mna_profiles_binned_tb %>%
  filter(Name %in% mna_samples) %>%
  group_by(start) %>%
  summarise(PercentSamples = 100*mean(isAmp),
            nSamples = sum(isAmp)) %>%
  ungroup() %>%
  ggplot(aes(x=start, y=PercentSamples)) +
  geom_line() + 
  xlim(15000000, 17000000) +
  theme_kons1() + 
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ylab("Samples [%]") + 
  xlab("chr2") + 
  ggtitle("Amplified Regions on chr2 (Depuydt et al. 2018)") +
  annotate(xmin=15307032, xmax=15701454, ymin=0, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # NBAS
  annotate(xmin=15731302, xmax=15771235, ymin=0, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # DDX1
  annotate(xmin=16060521, xmax=16076139, ymin=0, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # MYCNUT
  annotate(xmin=16080683, xmax=16087129, ymin=0, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # MYCN
  annotate(xmin=16190549, xmax=16225923, ymin=0, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # GACAT3
  annotate(xmin=16730727, xmax=16847599, ymin=0, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # FAM49A
  ggsave("~/Desktop/copynumber_HR_NB-master/AmplifiedOnChr2_MYCNFocus_All_GenesMarked.pdf",
         height = 3, width = 5, useDingbats=F)


# ------------------------------------------------------------------------------
# randomize with respect to MYCN
# ------------------------------------------------------------------------------

set.seed(42)
roi = GRanges(seqnames = "chr2",
              ranges = IRanges(start=16080683, end=16087129))
roi_length = width(roi)
randomized_profiles = list()
reps = 100
nrandomsamples = reps*240
for (i in 1:nrandomsamples){ # took about 10min
  random_patient = mna_profiles[[sample(names(mna_profiles), 1)]]
  random_patient = random_patient[width(random_patient)>roi_length]
  random_interval = random_patient[sample(length(random_patient), 1, prob=width(random_patient)-roi_length)]
  random_offset = sample.int(width(random_interval)-roi_length,1)
  absolute_offset = start(roi) - (start(random_interval) + random_offset)
  random_patient = trim(shift(random_patient, absolute_offset))
  randomized_profiles[[i]] = random_patient
}

library(parallel)
cores = detectCores()
# randomized_profiles = mclapply(1:nrandomsamples,
#                                function (i) {
#                                  set.seed(i)
#                                  random_patient = mna_profiles[[sample(names(mna_profiles), 1)]]
#                                  random_patient = random_patient[width(random_patient)>roi_length]
#                                  random_interval = random_patient[sample(length(random_patient), 1, prob=width(random_patient)-roi_length)]
#                                  random_offset = sample.int(width(random_interval)-roi_length,1)
#                                  absolute_offset = start(roi) - (start(random_interval) + random_offset)
#                                  random_patient = trim(shift(random_patient, absolute_offset))
#                                  return(random_patient)
#                                },
#                                mc.cores = cores)

randomized_profiles_cov = lapply(randomized_profiles, coverage)
randomized_profiles_binned = lapply(randomized_profiles_cov, function (this_profile_cov) binnedAverage(bins_focus, numvar = this_profile_cov, varname = "isAmp"))
#randomized_profiles_binned_tb = do.call(rbind, lapply(randomized_profiles_binned, as_tibble))
randomized_profiles_binned_tb = as_tibble(do.call(c, randomized_profiles_binned))

randomized_profiles_binned_tb %>%
  group_by(start) %>%
  summarise(PercentSamples = 100*mean(isAmp),
            nSamples = sum(isAmp)) %>%
  ungroup() %>%
  ggplot(aes(x=start, y=PercentSamples)) +
  geom_line() + 
  xlim(15000000, 17000000) +
  theme_kons1() + 
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ylab("Samples [%]") + 
  xlab("chr2") + 
  ggtitle("Amplified Regions on chr2 (Depuydt et al. 2018)") +
  annotate(xmin=15307032, xmax=15701454, ymin=0, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # NBAS
  annotate(xmin=15731302, xmax=15771235, ymin=0, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # DDX1
  annotate(xmin=16060521, xmax=16076139, ymin=0, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # MYCNUT
  annotate(xmin=16080683, xmax=16087129, ymin=0, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # MYCN
  annotate(xmin=16190549, xmax=16225923, ymin=0, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # GACAT3
  annotate(xmin=16730727, xmax=16847599, ymin=0, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) #+ # FAM49A
#ggsave("~/Desktop/copynumber_HR_NB-master/AmplifiedOnChr2__Randomized_MYCNFocus_All_GenesMarked.pdf",
#       height = 3, width = 5, useDingbats=F)


# Plot Real and Randomized Data in one plot
mna_profiles_binned_tb$RealOrRandomized = "Real"
randomized_profiles_binned_tb$RealOrRandomized = "Randomized"
bind_rows(mna_profiles_binned_tb, randomized_profiles_binned_tb) %>%
  group_by(start, RealOrRandomized) %>%
  summarise(PercentSamples = 100*mean(isAmp),
            nSamples = sum(isAmp)) %>%
  ungroup() %>%
  ggplot(aes(x=start, y=PercentSamples, linetype = RealOrRandomized)) +
  geom_line() + 
  xlim(15000000, 17000000) +
  theme_kons1() + 
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ylab("Samples [%]") + 
  xlab("chr2") + 
  scale_linetype_manual(values = c("Real"="solid", "Randomized"="dotted")) +
  ggtitle("Amplified Regions on chr2 (Depuydt et al. 2018)") +
  guides(linetype = F) +
  annotate(xmin=15307032, xmax=15701454, ymin=0, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # NBAS
  annotate(xmin=15731302, xmax=15771235, ymin=0, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # DDX1
  annotate(xmin=16060521, xmax=16076139, ymin=0, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # MYCNUT
  annotate(xmin=16080683, xmax=16087129, ymin=0, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # MYCN
  annotate(xmin=16190549, xmax=16225923, ymin=0, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # GACAT3
  annotate(xmin=16730727, xmax=16847599, ymin=0, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # FAM49A
  ggsave("/Volumes/Elements/MYCNAmplicon/Figures/AmplifiedOnChr2_RandomizedVsReal_MYCNFocus_All_GenesMarked.pdf",
         height = 3, width = 5, useDingbats=F)

# ------------------------------------------------------------------------------
# Calculate FC over chance
# ------------------------------------------------------------------------------

# Calcuale FC Over Chance
FCOverChance = bind_rows(mna_profiles_binned_tb, randomized_profiles_binned_tb) %>%
  group_by(start, RealOrRandomized) %>%
  summarise(PercentSamples = 100*mean(isAmp),
            nSamples = sum(isAmp)) %>%
  ungroup() %>%
  dplyr::select(-nSamples) %>% 
  tidyr::spread(RealOrRandomized, PercentSamples) %>%
  mutate(FC = (Real + 1)/(Randomized + 1))

# Create GenomicRanges object
FCOverChanceGR = FCOverChance %>%
  mutate(seqnames = "chr2", end = start + 9999, score=FC) %>%
  dplyr::select(seqnames, start, end, score) %>% 
  filter(!is.na(score)) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T) 
seqlevels(FCOverChanceGR) = seqlevels(BSgenome.Hsapiens.UCSC.hg19)
seqlengths(FCOverChanceGR) = seqlengths(BSgenome.Hsapiens.UCSC.hg19)

# Write bigwig file
FCOverChanceGR %>% 
  write_bigwig(file = "/Volumes/Elements/MYCNAmplicon/Figures/CopyNumber_FCOverChance.bigwig")

# Plot FC over Chance for MYCN locus
FCOverChance %>%
  filter(start > 15000000, start<17000000) %>% 
  ggplot(aes(x=start, y=FC)) +
  xlim(15000000, 17000000) + ylim(0.6,1.4) +
  geom_line() + 
  geom_hline(yintercept = 1, linetype="dashed") +
  theme_kons1() + 
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ylab("Fold Change Over Chance") + 
  xlab("chr2") + 
  annotate(xmin=15307032, xmax=15701454, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # NBAS
  annotate(xmin=15731302, xmax=15771235, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # DDX1
  annotate(xmin=16060521, xmax=16076139, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # MYCNUT
  annotate(xmin=16080683, xmax=16087129, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # MYCN
  annotate(xmin=16190549, xmax=16225923, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # GACAT3
  annotate(xmin=16730727, xmax=16847599, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # FAM49A
  ggsave("/Volumes/Elements/MYCNAmplicon/Figures/AmplifiedOnChr2_RandomizedVsReal_FCOverChance_MYCNFocus_All_GenesMarked.pdf",
         height = 3, width = 5, useDingbats=F)

# ------------------------------------------------------------------------------
# Make Fig 2a
# ------------------------------------------------------------------------------
genes_df = 
  data.frame(
    "start" = c(15307032, 15731302, 16060521, 16080683, 16190549, 16730727),
    "end" = c(15701454, 15771235, 16076139, 16087129, 16225923, 16847599),
    "name" = c("NBAS", "DDX1", "MYCNUT", "MYCN", "GACAT3", "FAM49A"),
    "class" = rep("gene", 6)
  )
genes.fig = genes_df %>%
  ggplot(aes(x=start, y=1)) +
  geom_rect(xmin=genes_df$start, xmax=genes_df$end, ymin=-Inf, ymax=Inf, color=NA, fill="black") + 
  theme_kons2() + 
  theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank()) +
  xlim(gene_of_interest_start-500000,gene_of_interest_end + 500000)
cRE = read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_nMNA_MYCN_Enhancers.bed") %>% as_tibble()
cRE.fig = cRE %>%
  ggplot(aes(x=start, y=score)) +
  geom_rect(aes(xmin=start, xmax=end), ymin=-Inf, ymax=Inf, color="firebrick3", fill="firebrick3") + 
  theme_kons2() + 
  theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank()) +
  xlim(gene_of_interest_start-500000,gene_of_interest_end + 500000)
cRE.fig 
cRE$class = "cRE"
cRE_genes_plot_df = 
  bind_rows(cRE, genes_df) %>% 
  mutate(dummy=1) %>% 
  dplyr::select(start, end, class, dummy) %>%
  mutate(this_color = ifelse(class == "gene", NA, "firebrick3"),
         this_fill = ifelse(class == "gene", "black", "firebrick3"))
cRE_genes.fig = 
  cRE_genes_plot_df %>%
  ggplot(aes(x=start, y=dummy)) +
  geom_rect(xmin=cRE_genes_plot_df$start, xmax=cRE_genes_plot_df$end, ymin=-Inf, ymax=Inf, color=cRE_genes_plot_df$this_color, fill=cRE_genes_plot_df$this_fill) + 
  theme_kons2() +
  theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank()) +
  xlim(gene_of_interest_start-500000,gene_of_interest_end + 500000)

real_profiles.fig =
  mna_profiles_binned_tb %>%
  filter(Name %in% mna_samples) %>%
  group_by(start) %>%
  summarise(PercentSamples = 100*mean(isAmp),
            nSamples = sum(isAmp)) %>%
  ungroup() %>%
  ggplot(aes(x=start, y=PercentSamples)) +
  geom_line() + 
  ylab("Samples [%]") + 
  xlim(gene_of_interest_start-500000,gene_of_interest_end + 500000) +
  theme_kons2() 

randomized_profiles.fig =
  randomized_profiles_binned_tb %>%
  group_by(start) %>%
  summarise(PercentSamples = 100*mean(isAmp),
            nSamples = sum(isAmp)) %>%
  ungroup() %>%
  ggplot(aes(x=start, y=PercentSamples)) +
  geom_line(linetype = "dashed") + 
  ylab("Samples [%]") + 
  xlim(gene_of_interest_start-500000,gene_of_interest_end + 500000) +
  theme_kons2() 

mna_profiles_binned_tb$RealOrRandomized = "Real"
randomized_profiles_binned_tb$RealOrRandomized = "Randomized"
real_and_randomized_profiles.fig = bind_rows(mna_profiles_binned_tb, randomized_profiles_binned_tb) %>%
  group_by(start, RealOrRandomized) %>%
  summarise(PercentSamples = 100*mean(isAmp),
            nSamples = sum(isAmp)) %>%
  ungroup() %>%
  ggplot(aes(x=start, y=PercentSamples, linetype = RealOrRandomized)) +
  geom_line() + 
  theme_kons1() + 
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ylab("Samples [%]") + 
  scale_linetype_manual(values = c("Real"="solid", "Randomized"="dotted")) +
  guides(linetype = F) +
  scale_y_continuous(breaks = c(40, 60, 80, 100), limits=c(40,100), expand=c(0,0)) +
  xlim(gene_of_interest_start-500000,gene_of_interest_end + 500000) +
  xlab("chr2") +
  theme_kons2() 

ggsave("/Volumes/Elements/MYCNAmplicon/Results/Fig2a_Nov23.pdf",
       egg::ggarrange(
         genes.fig+
           theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
         cRE.fig+
           theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
         real_and_randomized_profiles.fig,
         nrow = 3, 
         heights = c(0.05,0.05,1)),
       height=2, width=3, onefile = FALSE, useDingbats = F)



# ------------------------------------------------------------------------------
# calculate stats and FC for specific regions
# ------------------------------------------------------------------------------

# Calculate how many amplicons overlap with ROIs and how many would overlap by chance
# Compare Real vs. Random
enhancers = read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_nMNA_MYCN_Enhancers.bed") 
mna_profiles_rois = lapply(mna_profiles_cov, function (this_profile_cov) binnedAverage(enhancers, numvar = this_profile_cov, varname = "isAmp"))
mna_profiles_rois_tb = do.call(rbind, lapply(mna_profiles_rois, as_tibble))
mna_profiles_rois_tb$isAmp = as.numeric(mna_profiles_rois_tb$isAmp > 0) # ensure that partial overlap is also counted!

randomized_profiles_rois = lapply(randomized_profiles_cov, function (this_profile_cov) binnedAverage(enhancers, numvar = this_profile_cov, varname = "isAmp"))
randomized_profiles_rois_tb = do.call(rbind, lapply(randomized_profiles_rois, as_tibble))
randomized_profiles_rois_tb$isAmp = as.numeric(randomized_profiles_rois_tb$isAmp > 0) # ensure that partial overlap is also counted!
mna_profiles_rois_tb$RealOrRandomized = "Real"
randomized_profiles_rois_tb$RealOrRandomized = "Randomized"
bind_rows(mna_profiles_rois_tb, randomized_profiles_rois_tb) %>%
  group_by(name, RealOrRandomized) %>%
  summarise(PercentOfAmpliconsThatOverlap = 100*mean(isAmp),
            nSamples = sum(isAmp))
bind_rows(mna_profiles_rois_tb, randomized_profiles_rois_tb) %>%
  group_by(name, RealOrRandomized) %>%
  summarise(PercentOfAmpliconsThatOverlap = 100*mean(isAmp)) %>% 
  tidyr::spread(RealOrRandomized, PercentOfAmpliconsThatOverlap) %>%
  mutate(FC = (Real + 1)/(Randomized + 1))

randomization_df = list()
for (ri in 1:reps){
  startidx = 240*(ri-1) + 1
  endidx = 240*ri
  thisrep_randomized_profiles_cov = randomized_profiles_cov[startidx:endidx]
  randomized_profiles_rois = lapply(thisrep_randomized_profiles_cov, function (this_profile_cov) binnedAverage(enhancers, numvar = this_profile_cov, varname = "isAmp"))
  randomized_profiles_rois_tb = as_tibble(do.call(c, randomized_profiles_rois))
  randomized_profiles_rois_tb$isAmp = as.numeric(randomized_profiles_rois_tb$isAmp > 0) # ensure that partial overlap is also counted!
  randomization_df[[ri]] = randomized_profiles_rois_tb %>% group_by(name) %>% summarise(isAmp = mean(isAmp, na.rm=T)) %>% ungroup() %>% mutate(randomizationIdx = ri)
}
randomization_df = do.call(rbind, (randomization_df))

# get p value for e4

get_pval = function(ename){
  e4_real_value = mna_profiles_rois_tb %>% filter(name == ename) %>% .$isAmp %>% mean()
  e4_simulated_overlap = randomization_df %>% filter(name == ename) %>% .$isAmp
  e4_FC = e4_real_value / mean(e4_simulated_overlap, na.rm=T)
  e4_pval = (sum(e4_real_value <= e4_simulated_overlap, na.rm=T)+1) / (sum(!is.na(e4_simulated_overlap))+1)
  return(e4_pval)
}
get_random = function(ename){
  e4_simulated_overlap = randomization_df %>% filter(name == ename) %>% .$isAmp
  return(mean(e4_simulated_overlap, na.rm=T))
}
get_FC = function(ename){
  e4_real_value = mna_profiles_rois_tb %>% filter(name == ename) %>% .$isAmp %>% mean()
  e4_simulated_overlap = randomization_df %>% filter(name == ename) %>% .$isAmp
  e4_FC = e4_real_value / mean(e4_simulated_overlap, na.rm=T)
  e4_pval = (sum(e4_real_value <= e4_simulated_overlap, na.rm=T)+1) / (sum(!is.na(e4_simulated_overlap))+1)
  return(e4_FC)
}

mna_profiles_rois_tb_summary = mna_profiles_rois_tb %>%
  group_by(name) %>%
  summarise(isAmp = mean(isAmp, na.rm=T)) %>%
  ungroup() 
mna_profiles_rois_tb_summary$isRandomAmp= sapply(mna_profiles_rois_tb_summary$name, get_random)
mna_profiles_rois_tb_summary$pval = sapply(mna_profiles_rois_tb_summary$name, get_pval)
mna_profiles_rois_tb_summary$FC = sapply(mna_profiles_rois_tb_summary$name, get_FC)
mna_profiles_rois_tb_summary$pvalBH = p.adjust(mna_profiles_rois_tb_summary$pval, method = "BH")

# ------------------------------------------------------------------------------
# Upset plot for enhancer overlap (Fig. 2b)
# ------------------------------------------------------------------------------

library(UpSetR)

rois = read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_nMNA_MYCN_Enhancers.bed")
rois = rois[rois$name != "MYCNp"]

upset_list = list()
for (roi_idx in 1:length(rois)){
  this_vector = vector()
  for (sample_idx in 1:length(mna_profiles)){
    if (length(findOverlaps(rois[roi_idx], mna_profiles[sample_idx])) > 0) this_vector = c(this_vector, sample_idx)
  }
  upset_list[[rois[roi_idx]$name]] = this_vector
}

pdf(paste0("/Volumes/Elements/MYCNAmplicon/Results/EnhancersVsCopyNumber_UpsetPlot.pdf"), height=3, width=4, onefile = F)
upset(fromList(upset_list), 
      order.by="freq", 
      text.scale = 1, 
      point.size = 1,
      line.size = 0.25,
      mb.ratio = c(0.66, 0.33),
      mainbar.y.label = "Samples",
      sets.x.label = "",
      set_size.show = T)
dev.off()

# ------------------------------------------------------------------------------
# Upset plot for randomized overlap
# ------------------------------------------------------------------------------

library(UpSetR)

rois = read_bed("/Volumes/Elements/MYCNAmplicon/Figures/Boeva_nMNA_MYCN_Enhancers.bed")
rois = rois[rois$name != "MYCNp"]

randomized_profiles = GRangesList(randomized_profiles)
upset_list_randomized = list()
for (roi_idx in 1:length(rois)){
  this_vector = vector()
  for (sample_idx in 1:length(randomized_profiles)){
    if (length(findOverlaps(rois[roi_idx], randomized_profiles[sample_idx])) > 0) this_vector = c(this_vector, sample_idx)
  }
  upset_list_randomized[[rois[roi_idx]$name]] = this_vector
}

pdf(paste0("/Volumes/Elements/MYCNAmplicon/Figures/ROI_Randomized_upset.pdf"), height=4, width=4, onefile = F)
upset(fromList(upset_list_randomized), 
      order.by="freq", 
      text.scale = 1, 
      mainbar.y.label = "Samples",
      sets.x.label = "")
dev.off()

# ------------------------------------------------------------------------------
# Investigate those that do not have the enhancer
# ------------------------------------------------------------------------------
rois = read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_nMNA_MYCN_Enhancers.bed")
rois = rois[rois$name != "MYCNp"]

enhancer_overlap_df = data.frame(
  Sample = names(mna_profiles),
  nAmpliconParts = unname(sapply(mna_profiles, length))
)
for (roi_idx in 1:length(rois)){
  enhancer_overlap_df[,rois[roi_idx]$name] = NA
  for (sample_idx in 1:length(mna_profiles)){
    if (length(findOverlaps(rois[roi_idx], mna_profiles[sample_idx])) > 0){
      enhancer_overlap_df[sample_idx, rois[roi_idx]$name] = TRUE
    } else {
      enhancer_overlap_df[sample_idx, rois[roi_idx]$name] = FALSE
    }
  }
}

enhancer_overlap_df %>% 
  ggplot(aes(x=e4, y=nAmpliconParts)) +
  geom_jitter() + 
  theme_kons1()

# Do samples without the Enhancer have more amplicon parts?

enhancer_overlap_df %>% 
  group_by(e4, nAmpliconParts) %>% 
  summarise(n=dplyr::n()) %>%
  ggplot(aes(x=nAmpliconParts, y=n)) +
  geom_col()


n_e4 = sum(enhancer_overlap_df$e4)
n_non_e4 = sum(!enhancer_overlap_df$e4)

enhancer_overlap_df %>% 
  mutate(weight = ifelse(e4, 1/n_e4, 1/n_non_e4)) %>%
  group_by(e4, nAmpliconParts) %>% 
  summarise(perc=100*sum(weight)) %>%
  ungroup() %>% 
  mutate(e4 = factor(e4, levels=c("TRUE", "FALSE"))) %>% 
  ggplot(aes(x=nAmpliconParts, y=perc, fill=e4, color=e4, group=e4)) + 
  geom_col(position = position_dodge2(preserve = "single"), size=0.1) +
  guides(color=F, fill=F) +
  theme_kons2() +
  ylab("Samples [%]") + 
  xlab("Number of Fragments") +
  #ylim(0,100) +
  scale_color_manual(values=c("TRUE"="steelblue", "FALSE"="firebrick3")) +
  scale_fill_manual(values=c("TRUE"="steelblue", "FALSE"="firebrick3")) +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/MNAProfiles_nFragments_Perc.pdf", 
         height=1.5, width=1.5, useDingbats=F)
  
enhancer_overlap_df %>% 
  group_by(e4, nAmpliconParts) %>% 
  summarise(n=n()) %>%
  ungroup() %>% 
  mutate(e4 = factor(e4, levels=c("TRUE", "FALSE"))) %>% 
  ggplot(aes(x=nAmpliconParts, y=n, fill=e4, color=e4, group=e4)) + 
  geom_col(position = position_dodge2(preserve = "single"), size=0.1) +
  guides(color=F, fill=F) +
  theme_kons2() +
  ylab("Samples") + 
  xlab("Number of Fragments") + 
  scale_color_manual(values=c("TRUE"="steelblue", "FALSE"="firebrick3")) +
  scale_fill_manual(values=c("TRUE"="steelblue", "FALSE"="firebrick3")) +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/MNAProfiles_nFragments_abs.pdf", 
         height=1.5, width=1.5, useDingbats=F)
  
enhancer_overlap_df %>%
  ggplot(aes(x=nAmpliconParts)) +
  geom_histogram(binwidth = 1)
  


  ggplot(aes(y=nAmpliconParts)) +
  geom_jitter() + 
  theme_kons1()

# t-test
enhancer_overlap_df %>%
  group_by(e4) %>%
  summarise(n=dplyr::n(),
            meanNAmpliconParts = mean(nAmpliconParts))
t.test(nAmpliconParts ~ e4, data=enhancer_overlap_df)

# Fisher Test (One part vs. more parts)
con_table_df = enhancer_overlap_df %>% 
  mutate(moreThanOneAmpliconPart = nAmpliconParts>1) %>% 
  group_by(e4) %>% 
  summarise(n = dplyr::n(),
            CountOneAmpliconPart = sum(!moreThanOneAmpliconPart),
            CountMoreThanOneAmpliconPart = sum(moreThanOneAmpliconPart),
            PercentMoreThanOneAmpliconPart = 100*mean(moreThanOneAmpliconPart)) 

con_table = con_table_df %>%
  dplyr::select(CountOneAmpliconPart, CountMoreThanOneAmpliconPart) %>%
  as.matrix()
rownames(con_table) = ifelse(con_table_df$e4, "IncludesE4", "DoesNotIncludeE4")
print(con_table)
fisher.test(con_table)

# Now plot the percentages...
enhancer_overlap_df %>% 
  mutate(moreThanOneAmpliconPart = nAmpliconParts>1) %>% 
  group_by(e4) %>% 
  summarise(n=dplyr::n(), 
            nMoreThanOneAmpliconPart = sum(moreThanOneAmpliconPart),
            PercentMoreThanOneAmpliconPart = mean(moreThanOneAmpliconPart))
enhancer_overlap_df %>% 
  mutate(moreThanOneAmpliconPart = nAmpliconParts>1) %>% 
  group_by(e4) %>% 
  summarise(n=dplyr::n(), 
            nMoreThanOneAmpliconPart = sum(moreThanOneAmpliconPart),
            PercentMoreThanOneAmpliconPart = mean(moreThanOneAmpliconPart)) %>% 
  ggplot(aes(x=e4, y=PercentMoreThanOneAmpliconPart)) +
  geom_col() 


# ------------------------------------------------------------------------------
# Focus on the samples without e4 and only one part
# ------------------------------------------------------------------------------

# Maybe they have amplifications on other chromosomes?
mna_profiles_allchr = profiles_df %>%
  full_join(samples, by="Name") %>%
  filter(MYCNStatus == 1) %>%
  mutate(Sample = Name) %>% 
  mutate(chromosome = paste0("chr", as.character(chromosome))) %>%
  mutate(chromosome = as.character(chromosome)) %>% 
  filter(annotation == "ampl") %>%
  makeGRangesListFromDataFrame(keep.extra.columns = T,
                               split.field = "Name",
                               seqnames.field = "chromosome",
                               start.field = "min",
                               end.field = "max", 
                               seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg19))
seqlevels(mna_profiles_allchr, pruning.mode = "coarse") = seqlevels(BSgenome.Hsapiens.UCSC.hg19)
arrayDataIsMNA = vector()
for (sample_idx in 1:length(mna_profiles_allchr)){
  if (length(findOverlaps(mycn_granges, mna_profiles_allchr[sample_idx])) > 0) arrayDataIsMNA = c(arrayDataIsMNA, sample_idx)
}
mna_profiles_allchr = mna_profiles_allchr[arrayDataIsMNA]

bins_allchr = tileGenome(seqlengths = seqlengths(BSgenome.Hsapiens.UCSC.hg19), 
                         tilewidth = 10000,
                         cut.last.tile.in.chrom = T)
mna_profiles_cov_allchr = lapply(mna_profiles_allchr, coverage)
mna_profiles_binned_allchr = lapply(mna_profiles_cov_allchr, function (this_profile_cov) binnedAverage(bins_allchr, numvar = this_profile_cov, varname = "isAmp"))
for (i in 1:length(mna_profiles_binned_allchr)){
  mna_profiles_binned_allchr[[i]]$Name = names(mna_profiles_binned_allchr)[i]
}
mna_profiles_binned_tb_allchr = do.call(rbind, lapply(mna_profiles_binned_allchr, as_tibble))

mna_profiles_binned_tb_allchr %>%
  group_by(start) %>%
  summarise(PercentSamples = 100*mean(isAmp),
            nSamples = sum(isAmp)) %>%
  ungroup() %>%
  ggplot(aes(x=start, y=PercentSamples)) +
  geom_line() + 
  xlim(15000000, 17000000) +
  theme_kons1() + 
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ylab("Samples [%]") + 
  xlab("chr2")


# ------------------------------------------------------------------------------
# Investigate those that do not have the enhancer ALL CHR
# ------------------------------------------------------------------------------
rois = read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_nMNA_MYCN_Enhancers.bed")
rois = rois[rois$name != "MYCNp"]

enhancer_overlap_df = data.frame(
  Sample = names(mna_profiles_allchr),
  nAmpliconParts = unname(sapply(mna_profiles_allchr, length))
)
for (roi_idx in 1:length(rois)){
  enhancer_overlap_df[,rois[roi_idx]$name] = NA
  for (sample_idx in 1:length(mna_profiles_allchr)){
    if (length(findOverlaps(rois[roi_idx], mna_profiles_allchr[sample_idx])) > 0){
      enhancer_overlap_df[sample_idx, rois[roi_idx]$name] = TRUE
    } else {
      enhancer_overlap_df[sample_idx, rois[roi_idx]$name] = FALSE
    }
  }
}

enhancer_overlap_df %>% 
  ggplot(aes(x=e4, y=nAmpliconParts)) +
  geom_jitter() + 
  theme_kons1()

# Do samples without the Enhancer have more amplicon parts?

# t-test
enhancer_overlap_df %>%
  group_by(e4) %>%
  summarise(n=dplyr::n(),
            meanNAmpliconParts = mean(nAmpliconParts))
t.test(nAmpliconParts ~ e4, data=enhancer_overlap_df)

# Fisher Test (One part vs. more parts)
con_table_df = enhancer_overlap_df %>% 
  mutate(moreThanOneAmpliconPart = nAmpliconParts>1) %>% 
  group_by(e4) %>% 
  summarise(n = dplyr::n(),
            CountOneAmpliconPart = sum(!moreThanOneAmpliconPart),
            CountMoreThanOneAmpliconPart = sum(moreThanOneAmpliconPart),
            PercentMoreThanOneAmpliconPart = 100*mean(moreThanOneAmpliconPart)) 

con_table = con_table_df %>%
  dplyr::select(CountOneAmpliconPart, CountMoreThanOneAmpliconPart) %>%
  as.matrix()
rownames(con_table) = ifelse(con_table_df$e4, "IncludesE4", "DoesNotIncludeE4")
sink(file="/Volumes/Elements/MYCNAmplicon/Results/E4NotIncluded_NFragments.txt")
print(con_table) 
fisher.test(con_table)
sink()

# Now plot the percentages...
enhancer_overlap_df %>% 
  mutate(moreThanOneAmpliconPart = nAmpliconParts>1) %>% 
  group_by(e4) %>% 
  summarise(n=dplyr::n(), 
            nMoreThanOneAmpliconPart = sum(moreThanOneAmpliconPart),
            PercentMoreThanOneAmpliconPart = mean(moreThanOneAmpliconPart))
enhancer_overlap_df %>% 
  mutate(moreThanOneAmpliconPart = nAmpliconParts>1) %>% 
  group_by(e4) %>% 
  summarise(n=dplyr::n(), 
            nMoreThanOneAmpliconPart = sum(moreThanOneAmpliconPart),
            PercentMoreThanOneAmpliconPart = mean(moreThanOneAmpliconPart)) %>% 
  ggplot(aes(x=e4, y=PercentMoreThanOneAmpliconPart)) +
  geom_col() 


# ------------------------------------------------------------------------------
# What about samples that do not include e4
# ------------------------------------------------------------------------------

e4less_samples = 
  enhancer_overlap_df %>% 
  filter(!e4) %>% 
  .$Sample

mna_profiles[e4less_samples] %>%
  write_bed("/Volumes/Elements/MYCNAmplicon/Results/MNAProfilesChr2_E4less.bed")

e4less_samples_mna_profiles = mna_profiles_allchr[e4less_samples] %>% unlist()
e4less_samples_mna_profiles$Sample = names(e4less_samples_mna_profiles)
e4less_samples_mna_profiles = unname(e4less_samples_mna_profiles)

e3less_samples = enhancer_overlap_df %>% 
  filter(!e3) %>% 
  .$Sample
e3less_samples_mna_profiles = mna_profiles_allchr[e3less_samples] %>% unlist()
e3less_samples_mna_profiles$Sample = names(e3less_samples_mna_profiles)
e3less_samples_mna_profiles = unname(e3less_samples_mna_profiles)

enhancer_overlap_df %>% dplyr::select(Sample, e3, e4) %>% View

e3ANDe4less_samples = enhancer_overlap_df %>% 
  filter(!e3 & !e4) %>% 
  .$Sample
e3ANDe4less_samples_mna_profiles = mna_profiles_allchr[e3ANDe4less_samples] %>% unlist()
e3ANDe4less_samples_mna_profiles$Sample = names(e3ANDe4less_samples_mna_profiles)
e3ANDe4less_samples_mna_profiles = unname(e3ANDe4less_samples_mna_profiles)

e3ORe4less_samples = enhancer_overlap_df %>% 
  filter(!e3 | !e4) %>% 
  .$Sample
e3ORe4less_samples_mna_profiles = mna_profiles_allchr[e3ORe4less_samples] %>% unlist()
e3ORe4less_samples_mna_profiles$Sample = names(e3ORe4less_samples_mna_profiles)
e3ORe4less_samples_mna_profiles = unname(e3ORe4less_samples_mna_profiles)

# which_enhancers_how_often = CRC_enhancers %>%
#   group_by_overlaps(e4less_samples_mna_profiles) %>%
#   plyranges::summarise(n=plyranges::n()) %>% 
#   as_tibble
# CRC_enhancers$howoftencoamplified = NA
# CRC_enhancers[which_enhancers_how_often$query]$howoftencoamplified = which_enhancers_how_often$n
# CRC_enhancers %>% as_tibble() %>% View
# CRC_enhancers = read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_lowMYCNExpr_H3K27ac_AggregatePeaks_CRCFactorPositiveOnly.bed")
# CRC_enhancers %>% filter(seqnames == "chr2")

e4less_mna_profiles = mna_profiles[e4less_samples]
save(e4less_mna_profiles, file="/Volumes/Elements/MYCNAmplicon/Results/e4less_MNAProfiles.Rdata")

e3less_mna_profiles = mna_profiles[e3less_samples]
save(e3less_mna_profiles, file="/Volumes/Elements/MYCNAmplicon/Results/e3less_MNAProfiles.Rdata")

e3ORe4less_mna_profiles = mna_profiles[e3ORe4less_samples]
save(e3ORe4less_mna_profiles, file="/Volumes/Elements/MYCNAmplicon/Results/e3ORe4less_MNAProfiles.Rdata")

e3ANDe4less_mna_profiles = mna_profiles[e3ANDe4less_samples]
save(e3ANDe4less_mna_profiles, file="/Volumes/Elements/MYCNAmplicon/Results/e3ANDe4less_MNAProfiles.Rdata")

# ------------------------------------------------------------------------------
# super enhancer overlap
# ------------------------------------------------------------------------------

lowmycn_crc_se = read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_lowMYCN_CRCdriven_SE.bed")
lowmycn_se = read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_lowMYCN_SE.bed")
nonmna_se = read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_noMNA_SE.bed")

no_lowMYCN_SE = vector()
no_nMNA_SE = vector()
no_lowMYCN_crc_SE = vector()
for (sample_idx in 1:length(mna_profiles_allchr)){
  if (length(findOverlaps(lowmycn_crc_se, mna_profiles_allchr[sample_idx])) == 0) no_lowMYCN_crc_SE = c(no_lowMYCN_crc_SE, sample_idx)
  if (length(findOverlaps(lowmycn_se, mna_profiles_allchr[sample_idx])) == 0) no_lowMYCN_SE = c(no_lowMYCN_SE, sample_idx)
  if (length(findOverlaps(nonmna_se, mna_profiles_allchr[sample_idx])) == 0) no_nMNA_SE = c(no_nMNA_SE, sample_idx)
}

no_lowMYCN_crc_SE = names(mna_profiles_allchr)[no_lowMYCN_crc_SE]
length(no_lowMYCN_crc_SE)

no_lowMYCN_SE = names(mna_profiles_allchr)[no_lowMYCN_SE]
length(no_lowMYCN_SE)

no_nMNA_SE = names(mna_profiles_allchr)[no_nMNA_SE]
length(no_nMNA_SE)

# e4less samples
e4lesssamples_no_lowMYCN_SE = vector()
e4lesssamples_no_nMNA_SE = vector()
e4lesssamples_no_lowMYCN_crc_SE = vector()
for (sample in unique(e4less_samples_mna_profiles$Sample)){
  if (length(findOverlaps(lowmycn_crc_se, e4less_samples_mna_profiles[e4less_samples_mna_profiles$Sample == sample])) == 0) e4lesssamples_no_lowMYCN_crc_SE = c(e4lesssamples_no_lowMYCN_crc_SE, sample)
  if (length(findOverlaps(lowmycn_se, e4less_samples_mna_profiles[e4less_samples_mna_profiles$Sample == sample])) == 0) e4lesssamples_no_lowMYCN_SE = c(e4lesssamples_no_lowMYCN_SE, sample)
  if (length(findOverlaps(nonmna_se, e4less_samples_mna_profiles[e4less_samples_mna_profiles$Sample == sample])) == 0) e4lesssamples_no_nMNA_SE = c(e4lesssamples_no_nMNA_SE, sample)
}
length(e4lesssamples_no_lowMYCN_crc_SE)
length(e4lesssamples_no_lowMYCN_SE)
length(e4lesssamples_no_nMNA_SE)

# e3less samples
e3lesssamples_no_lowMYCN_SE = vector()
e3lesssamples_no_nMNA_SE = vector()
e3lesssamples_no_lowMYCN_crc_SE = vector()
for (sample in unique(e3less_samples_mna_profiles$Sample)){
  if (length(findOverlaps(lowmycn_crc_se, e3less_samples_mna_profiles[e3less_samples_mna_profiles$Sample == sample])) == 0) e3lesssamples_no_lowMYCN_crc_SE = c(e3lesssamples_no_lowMYCN_crc_SE, sample)
  if (length(findOverlaps(lowmycn_se, e3less_samples_mna_profiles[e3less_samples_mna_profiles$Sample == sample])) == 0) e3lesssamples_no_lowMYCN_SE = c(e3lesssamples_no_lowMYCN_SE, sample)
  if (length(findOverlaps(nonmna_se, e3less_samples_mna_profiles[e3less_samples_mna_profiles$Sample == sample])) == 0) e3lesssamples_no_nMNA_SE = c(e3lesssamples_no_nMNA_SE, sample)
}
length(e3lesssamples_no_lowMYCN_crc_SE)
length(e3lesssamples_no_lowMYCN_SE)
length(e3lesssamples_no_nMNA_SE)

# e3ANDe4less samples
e3ANDe4lesssamples_no_lowMYCN_SE = vector()
e3ANDe4lesssamples_no_nMNA_SE = vector()
e3ANDe4lesssamples_no_lowMYCN_crc_SE = vector()
for (sample in unique(e3ANDe4less_samples_mna_profiles$Sample)){
  if (length(findOverlaps(lowmycn_crc_se, e3ANDe4less_samples_mna_profiles[e3ANDe4less_samples_mna_profiles$Sample == sample])) == 0) e3ANDe4lesssamples_no_lowMYCN_crc_SE = c(e3ANDe4lesssamples_no_lowMYCN_crc_SE, sample)
  if (length(findOverlaps(lowmycn_se, e3ANDe4less_samples_mna_profiles[e3ANDe4less_samples_mna_profiles$Sample == sample])) == 0) e3ANDe4lesssamples_no_lowMYCN_SE = c(e3ANDe4lesssamples_no_lowMYCN_SE, sample)
  if (length(findOverlaps(nonmna_se, e3ANDe4less_samples_mna_profiles[e3ANDe4less_samples_mna_profiles$Sample == sample])) == 0) e3ANDe4lesssamples_no_nMNA_SE = c(e3ANDe4lesssamples_no_nMNA_SE, sample)
}
length(e3ANDe4lesssamples_no_lowMYCN_crc_SE)
length(e3ANDe4lesssamples_no_lowMYCN_SE)
length(e3ANDe4lesssamples_no_nMNA_SE)

# ------------------------------------------------------------------------------
# e4 less samples
# ------------------------------------------------------------------------------
library(regioneR)

blacklist_fname = "/Volumes/Elements/MYCNAmplicon/Data/hg19-blacklist.v2.bed"
encode_mask = data.table::fread(blacklist_fname) %>%
  as_tibble() %>%
  mutate(seqnames = V1, start = V2, end=V3) %>%
  makeGRangesFromDataFrame()

e4less_samples_mna_profiles_withoutMYCNlocus = 
  e4less_samples_mna_profiles %>%
  filter(seqnames == "chr2")

lowmycn_crc_se = lowmycn_crc_se %>%
  filter(start != 16177799, start != 16355899)

e4less_samples_mna_profiles_withoutMYCNlocus$doesOverlap = FALSE
doesoverlap = unique(queryHits(findOverlaps(e4less_samples_mna_profiles_withoutMYCNlocus, lowmycn_crc_se)))
e4less_samples_mna_profiles_withoutMYCNlocus[doesoverlap]$doesOverlap = TRUE
real_percent_overlap = 
  e4less_samples_mna_profiles_withoutMYCNlocus %>%
  as_tibble() %>%
  group_by(Sample) %>%
  summarise(anyoverlap = any(doesOverlap)) %>%
  ungroup() %>%
  summarise(percentHits = mean(anyoverlap, na.rm=T)) %>%
  .$percentHits
print(real_percent_overlap)

percent_hits = list()
for (i in 1:100){
  set.seed(i)
  e4less_samples_mna_profiles_randomized = 
    randomizeRegions(
      e4less_samples_mna_profiles_withoutMYCNlocus, 
      per.chromosome=T,
      genome="hg19")
  e4less_samples_mna_profiles_randomized$Sample = e4less_samples_mna_profiles_withoutMYCNlocus$Sample
  e4less_samples_mna_profiles_randomized$doesOverlap = FALSE
  doesoverlap = unique(queryHits(findOverlaps(e4less_samples_mna_profiles_randomized, nonmna_se)))
  e4less_samples_mna_profiles_randomized[doesoverlap]$doesOverlap = TRUE
  
  percent_hits[[i]] = 
    e4less_samples_mna_profiles_randomized %>%
    as_tibble() %>%
    group_by(Sample) %>%
    summarise(anyoverlap = any(doesOverlap)) %>%
    ungroup() %>%
    summarise(percentHits = mean(anyoverlap, na.rm=T))
}
percent_hits = unlist(percent_hits)
percent_hits = unname(percent_hits)
print(percent_hits)

pval = (sum(percent_hits>=real_percent_overlap)+1) / (length(percent_hits) + 1)
print(pval)


# ------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------  

e3ANDe4less_samples_mna_profiles_withoutMYCNlocus = 
  e3ANDe4less_samples_mna_profiles %>%
  filter(seqnames=="chr2")

e3ANDe4less_samples_mna_profiles_withoutMYCNlocus = mna_profiles_allchr %>% unlist()
e3ANDe4less_samples_mna_profiles_withoutMYCNlocus$Sample = names(e3ANDe4less_samples_mna_profiles_withoutMYCNlocus)
e3ANDe4less_samples_mna_profiles_withoutMYCNlocus = unname(e3ANDe4less_samples_mna_profiles_withoutMYCNlocus)

lowmycn_crc_se = lowmycn_crc_se %>%
  filter(start != 16177799, start != 16355899)

nonmna_se = nonmna_se %>%
  filter(start != 16177799, start != 16355899)

e3ANDe4less_samples_mna_profiles_withoutMYCNlocus$doesOverlap = FALSE
doesoverlap = unique(queryHits(findOverlaps(e3ANDe4less_samples_mna_profiles_withoutMYCNlocus, lowmycn_crc_se, minoverlap=10000)))
e3ANDe4less_samples_mna_profiles_withoutMYCNlocus[doesoverlap]$doesOverlap = TRUE
real_percent_overlap = 
  e3ANDe4less_samples_mna_profiles_withoutMYCNlocus %>%
  as_tibble() %>%
  group_by(Sample) %>%
  summarise(anyoverlap = any(doesOverlap)) %>%
  ungroup() %>%
  summarise(percentHits = mean(anyoverlap, na.rm=T)) %>%
  .$percentHits
print(real_percent_overlap)

percent_hits = list()
for (i in 1:5){
  set.seed(i)
  e3ANDe4less_samples_mna_profiles_randomized = 
    randomizeRegions(e3ANDe4less_samples_mna_profiles_withoutMYCNlocus, 
                     per.chromosome=T,
                     genome="hg19",
                     mask = encode_mask)
  e3ANDe4less_samples_mna_profiles_randomized$Sample = e3ANDe4less_samples_mna_profiles_withoutMYCNlocus$Sample
  e3ANDe4less_samples_mna_profiles_randomized$doesOverlap = FALSE
  doesoverlap = unique(queryHits(findOverlaps(e3ANDe4less_samples_mna_profiles_randomized, lowmycn_crc_se, minoverlap=10000)))
  e3ANDe4less_samples_mna_profiles_randomized[doesoverlap]$doesOverlap = TRUE
  
  percent_hits[[i]] = 
    e3ANDe4less_samples_mna_profiles_randomized %>%
    as_tibble() %>%
    group_by(Sample) %>%
    summarise(anyoverlap = any(doesOverlap)) %>%
    ungroup() %>%
    summarise(percentHits = mean(anyoverlap, na.rm=T))
}
percent_hits = unlist(percent_hits)
percent_hits = unname(percent_hits)
print(percent_hits)

pval = (sum(percent_hits>=real_percent_overlap)+1) / (length(percent_hits) + 1)
print(pval)


# ------------------------------------------------------------------------------

mna_profiles_allchr_df = mna_profiles_allchr %>% unlist()
mna_profiles_allchr_df$Sample = names(mna_profiles_allchr_df)
mna_profiles_allchr_df = unname(mna_profiles_allchr_df)
mna_profiles_allchr_df = as_tibble(mna_profiles_allchr_df)

mna_profiles_allchr_df %>%
  mutate(width = end-start) %>% 
  group_by(Sample) %>% 
  summarise(widthsum = sum(width, na.rm=T),
            nFragments = dplyr::n(),
            moreThanOneFragment = nFragments > 1) %>% 
  ungroup() %>%
  mutate(lacksE4 = Sample %in% e4less_samples) %>%
  ggplot(aes(x = lacksE4, y=nFragments)) + 
  geom_jitter()

# ------------------------------------------------------------------------------
# Clinical analysis
# ------------------------------------------------------------------------------

library(survival)
library(survminer)

# First only amplicon fragments
for (i in 1:length(mna_profiles)){
  mna_profiles[[i]]$Sample = names(mna_profiles)[i]
}
for (i in 1:length(mna_profiles_allchr)){
  mna_profiles_allchr[[i]]$Sample = names(mna_profiles_allchr)[i]
}

survival_data = 
  mna_profiles_allchr %>% 
  unlist() %>%
  as_tibble() %>% 
  group_by(Sample, OS, OStime) %>%
  summarise(nFragment = dplyr::n()) %>% 
  dplyr::select(Sample,OS,OStime,nFragment) %>%
  distinct() %>%
  mutate(MoreThanOneFragment = nFragment>1) %>%
  mutate(MoreThanOneFragment = factor(MoreThanOneFragment, levels=c("TRUE", "FALSE"))) 
survival_data$OSSurv = Surv(time=survival_data$OStime, event=survival_data$OS)
fit <- survfit(OSSurv ~ MoreThanOneFragment, data = survival_data)
pdf("/Volumes/Elements/MYCNAmplicon/Results/Clinical_MNA_moreThanOneAmplifiedFragment.pdf",
    height=4.5, width=5, useDingbats=F, onefile = F)
ggsurvplot(fit, data = survival_data, pval = TRUE, palette="Set1", risk.table = T, legend.title = "", ylab="Overall Survival", xlab = "Days",
           theme = 
             theme_survminer() + 
             theme(text = element_text(family="Helvetica", size=6),
                   title = element_text(family="Helvetica", size=6),
                   axis.text = element_text(family="Helvetica", size=6),
                   axis.title = element_text(family="Helvetica", size=6)
             ))
dev.off()
# 240 patients, 4 patients do not have survival data

survival_data %>%
  group_by(MoreThanOneFragment) %>% 
  summarise(n=n_distinct(Sample))
# 91 / (149+91)

survival_data = survival_data %>%
  mutate(MoreThanOneFragment = factor(MoreThanOneFragment, levels=c("FALSE", "TRUE"))) 
fit.coxph <- coxph(OSSurv ~ MoreThanOneFragment, data = survival_data)
pdf("/Volumes/Elements/MYCNAmplicon/Results/Clinical_MNA_moreThanOneAmplifiedFragment_Forest.pdf",
    height=4.5, width=5, useDingbats=F, onefile = F)
ggforest(fit.coxph, data = survival_data)
dev.off()

e4less_samples = 
  enhancer_overlap_df %>% 
  filter(!e4) %>% 
  .$Sample
survival_data = survival_data %>%
  mutate(LacksE4 = Sample %in% e4less_samples)
survival_data$OSSurv = Surv(time=survival_data$OStime, event=survival_data$OS)
fit <- survfit(OSSurv ~ LacksE4, data = survival_data)
pdf("/Volumes/Elements/MYCNAmplicon/Results/Clinical_MNA_LacksE4.pdf",
    height=4.5, width=5, useDingbats=F, onefile = F)
ggsurvplot(fit, data = survival_data, pval = TRUE, palette="Set1", risk.table = T, legend.title = "", ylab="Overall Survival", xlab = "Days",
           theme = 
             theme_survminer() + 
             theme(text = element_text(family="Helvetica", size=6),
                   title = element_text(family="Helvetica", size=6),
                   axis.text = element_text(family="Helvetica", size=6),
                   axis.title = element_text(family="Helvetica", size=6)
             ))
dev.off()

fit.coxph <- coxph(OSSurv ~ LacksE4, data = survival_data)
pdf("/Volumes/Elements/MYCNAmplicon/Results/Clinical_MNA_LacksE4_Forest.pdf",
    height=4.5, width=5, useDingbats=F, onefile = F)
ggforest(fit.coxph, data = survival_data)
dev.off()

odc1 = GRanges(
  seqnames = c(              "chr2"  ),
  ranges = IRanges(start = c(10580094),
                   end =   c(10588630))
)

distance_odc1_to_mycn = 16080683-10588630
odc1$GeneName = "ODC1"
survival_data = 
  mna_profiles_allchr %>% 
  unlist() %>%
  join_overlap_left(odc1, minoverlap = 10588630-10580094) %>%
  as_tibble() %>% 
  group_by(Sample, OS, OStime) %>%
  summarise(nFragment = dplyr::n(),
            ODC1CoAmpl = any(!is.na(GeneName))) %>% 
  mutate(ODC1CoAmpl = factor(ODC1CoAmpl, levels=c("TRUE", "FALSE"))) %>% 
  dplyr::select(Sample,OS,OStime,ODC1CoAmpl) %>%
  distinct()

survival_data %>%
  group_by(ODC1CoAmpl) %>% 
  summarise(n=n_distinct(Sample))
#21/(219+21)

survival_data$OSSurv = Surv(time=survival_data$OStime, event=survival_data$OS)
fit <- survfit(OSSurv ~ ODC1CoAmpl, data = survival_data)
pdf("/Volumes/Elements/MYCNAmplicon/Results/Clinical_MNA_ODC1CoAmpl.pdf",
    height=4.5, width=5, useDingbats=F, onefile = F)
ggsurvplot(fit, data = survival_data, pval = TRUE, palette="Set1", risk.table = T, legend.title = "", ylab="Overall Survival", xlab = "Days",
           theme = 
             theme_survminer() + 
             theme(text = element_text(family="Helvetica", size=6),
                   title = element_text(family="Helvetica", size=6),
                   axis.text = element_text(family="Helvetica", size=6),
                   axis.title = element_text(family="Helvetica", size=6)
             ))
dev.off()

# 240 patients, 4 patients do not have survival data
survival_data = survival_data %>%
  mutate(ODC1CoAmpl = factor(ODC1CoAmpl, levels=c("FALSE", "TRUE"))) 
fit.coxph <- coxph(OSSurv ~ ODC1CoAmpl, data = survival_data)
pdf("/Volumes/Elements/MYCNAmplicon/Results/Clinical_MNA_ODC1CoAmpl_Forest.pdf",
    height=4.5, width=5, useDingbats=F, onefile = F)
ggforest(fit.coxph, data = survival_data)
dev.off()
odc1coampl_samples = survival_data[survival_data$ODC1CoAmpl == "TRUE","Sample"] %>% .$Sample

greb1 = GRanges(
  seqnames = c(              "chr2"  ),
  ranges = IRanges(start = c(11674242),
                   end =   c(11782914))
)
greb1$GeneName = "GREB1"
survival_data = 
  mna_profiles_allchr %>% 
  unlist() %>%
  join_overlap_left(greb1, minoverlap = 11782914-11674242) %>%
  as_tibble() %>% 
  group_by(Sample, OS, OStime) %>%
  summarise(nFragment = dplyr::n(),
            GREB1CoAmpl = any(!is.na(GeneName))) %>% 
  dplyr::select(Sample,OS,OStime,GREB1CoAmpl) %>%
  mutate(GREB1CoAmpl = factor(GREB1CoAmpl, levels=c("TRUE", "FALSE"))) %>% 
  distinct()
survival_data %>%
  group_by(GREB1CoAmpl) %>% 
  summarise(n=n_distinct(Sample))
survival_data$OSSurv = Surv(time=survival_data$OStime, event=survival_data$OS)
fit <- survfit(OSSurv ~ GREB1CoAmpl, data = survival_data)
pdf("/Volumes/Elements/MYCNAmplicon/Results/Clinical_MNA_GREB1CoAmpl.pdf",
    height=4.5, width=5, useDingbats=F, onefile = F)
ggsurvplot(fit, data = survival_data, pval = TRUE, palette="Set1", risk.table = T, legend.title = "", ylab="Overall Survival", xlab = "Days",
           theme = 
             theme_survminer() + 
             theme(text = element_text(family="Helvetica", size=6),
                   title = element_text(family="Helvetica", size=6),
                   axis.text = element_text(family="Helvetica", size=6),
                   axis.title = element_text(family="Helvetica", size=6)
             ))
dev.off()

alk = GRanges(
  seqnames = c(              "chr2"  ),
  ranges = IRanges(start = c(29415640),
                   end =   c(30144477))
)
alk$GeneName = "ALK"
survival_data = 
  mna_profiles_allchr %>% 
  unlist() %>%
  join_overlap_left(alk, minoverlap = 30144477-29415640) %>%
  as_tibble() %>% 
  group_by(Sample, OS, OStime) %>%
  summarise(nFragment = dplyr::n(),
            ALKCoAmpl = any(!is.na(GeneName))) %>% 
  dplyr::select(Sample,OS,OStime,ALKCoAmpl) %>%
  mutate(ALKCoAmpl = factor(ALKCoAmpl, levels=c("TRUE", "FALSE"))) %>% 
  distinct()
survival_data %>%
  group_by(ALKCoAmpl) %>% 
  summarise(n=n_distinct(Sample))
survival_data$OSSurv = Surv(time=survival_data$OStime, event=survival_data$OS)
fit <- survfit(OSSurv ~ ALKCoAmpl, data = survival_data)
pdf("/Volumes/Elements/MYCNAmplicon/Results/Clinical_MNA_ALKCoAmpl.pdf",
    height=4.5, width=5, useDingbats=F, onefile = F)
ggsurvplot(fit, data = survival_data, pval = TRUE, palette="Set1", risk.table = T, legend.title = "", ylab="Overall Survival", xlab = "Days",
           theme = 
             theme_survminer() + 
             theme(text = element_text(family="Helvetica", size=6),
                   title = element_text(family="Helvetica", size=6),
                   axis.text = element_text(family="Helvetica", size=6),
                   axis.title = element_text(family="Helvetica", size=6)
             ))
dev.off()


survival_data = survival_data %>%
  mutate(ALKCoAmpl = factor(ALKCoAmpl, levels=c("FALSE", "TRUE"))) 
fit.coxph <- coxph(OSSurv ~ ALKCoAmpl, data = survival_data)
pdf("/Volumes/Elements/MYCNAmplicon/Results/Clinical_MNA_ALKCoAmpl_Forest.pdf",
    height=4.5, width=5, useDingbats=F, onefile = F)
ggforest(fit.coxph, data = survival_data)
dev.off()

odc1coampl_samples = survival_data[survival_data$ODC1CoAmpl == "TRUE","Sample"] %>% .$Sample
alkcoampl_samples = survival_data[survival_data$ALKCoAmpl == "TRUE","Sample"] %>% .$Sample

# ------------------------------------------------------------------------------
# Playground Clinical analysis
# ------------------------------------------------------------------------------

survival_data = survival_data %>%
  dplyr::select(Sample, OS, OStime) %>% 
  mutate(ALKCoAmpl = Sample %in% alkcoampl_samples, 
         ODC1CoAmpl = Sample %in% odc1coampl_samples) 
survival_data$OSSurv = Surv(time=survival_data$OStime, event=survival_data$OS)
fit.coxph <- coxph(OSSurv ~ ODC1CoAmpl + ALKCoAmpl, data = survival_data)
ggforest(fit.coxph, data = survival_data)

# mna_profiles_allchr_stitched = list()
# for (sample in names(mna_profiles_allchr)){
#   mna_profiles_allchr_stitched[[sample]] = reduce(mna_profiles_allchr[[sample]], min.gapwidth=100000) %>%
#     mutate(Sample = mna_profiles_allchr[[sample]]$Sample[1],
#            OS = mna_profiles_allchr[[sample]]$OS[1],
#            OStime = mna_profiles_allchr[[sample]]$OStime[1])
# 
# }
# mna_profiles_allchr_stitched = GRangesList(mna_profiles_allchr_stitched)

mna_profiles_allchr_stitched = list()
for (sample in names(mna_profiles_allchr)){
  mna_profiles_allchr_stitched[[sample]] = reduce(mna_profiles_allchr[[sample]], min.gapwidth=100000) %>%
    mutate(Sample = mna_profiles_allchr[[sample]]$Sample[1],
           OS = mna_profiles_allchr[[sample]]$OS[1],
           OStime = mna_profiles_allchr[[sample]]$OStime[1])
  
}
mna_profiles_allchr_stitched = GRangesList(mna_profiles_allchr_stitched)
survival_data = 
  mna_profiles_allchr %>% 
  unlist() %>%
  filter(seqnames=="chr2") %>% 
  as_tibble() %>% 
  group_by(Sample, OS, OStime) %>%
  summarise(nFragment = dplyr::n()) %>% 
  dplyr::select(Sample,OS,OStime,nFragment) %>%
  distinct() %>%
  mutate(FragmentDichot = nFragment>1)
survival_data$OSSurv = Surv(time=survival_data$OStime, event=survival_data$OS)
fit <- survfit(OSSurv ~ FragmentDichot, data = survival_data)
ggsurvplot(fit, data = survival_data, pval = TRUE, palette="Set1", risk.table = T, legend.title = "", ylab="Overall Survival", xlab = "Days")

goisss = GRanges(
  seqnames = c(              "chr2"  ),
  ranges = IRanges(start = c(29415640),
                   end =   c(30144477))
)
goisss$GeneName = "ALK"

goisss = GRanges(
  seqnames = c(              "chr2"  ),
  ranges = IRanges(start = c(10580094),
                   end =   c(10588630))
)
goisss$GeneName = "ODC1"


goisss = GRanges(
  seqnames = c(              "chr2"  ),
  ranges = IRanges(start = c(16190549),
                   end =   c(16225923))
)
goisss$GeneName = "GACAT3"


goisss = GRanges(
  seqnames = c(              "chr2"  ),
  ranges = IRanges(start = c(16375277),
                   end =   c(16385344))
)
goisss$GeneName = "e4"


goisss = GRanges(
  seqnames = c(              "chr2"  ),
  ranges = IRanges(start = c(24714783),
                   end =   c(24993571))
)
goisss$GeneName = "NCOA1"

goisss = GRanges(
  seqnames = c(              "chr2"  ),
  ranges = IRanges(start = c(6000000),
                   end =   c(7000000))
)
goisss$GeneName = "WeirdRegion"


goisss = GRanges(
  seqnames = c(              "chr12"  ),
  ranges = IRanges(start = c(58140616),
                   end =   c(58146228))
)
goisss$GeneName = "WeirdRegion"

goisss = GRanges(
  seqnames = c(              "chr12"  ),
  ranges = IRanges(start = c(69199952),
                   end =   c(69241324))
)
goisss$GeneName = "WeirdRegion"

goisss = GRanges(
  seqnames = c(              "chr2"  ),
  ranges = IRanges(start = c(6645149),
                   end =   c(6668251)))
goisss$GeneName = "SE"
#chr2:5832799 5841517

survival_data = 
  mna_profiles_allchr %>% 
  unlist() %>%
  join_overlap_left(goisss) %>%
  as_tibble() %>% 
  group_by(Sample, OS, OStime) %>%
  summarise(nFragment = dplyr::n(),
            doesHaveALK = any(!is.na(GeneName))) %>% 
  dplyr::select(Sample,OS,OStime,nFragment, doesHaveALK) %>%
  distinct() %>%
  mutate(lacksE4 = doesHaveALK)
survival_data$OSSurv = Surv(time=survival_data$OStime, event=survival_data$OS)
fit <- survfit(OSSurv ~ lacksE4, data = survival_data)
ggsurvplot(fit, data = survival_data, pval = TRUE, palette="Set1", risk.table = T, legend.title = "", ylab="Overall Survival", xlab = "Days")

survival_data = 
  mna_profiles_allchr %>% 
  unlist() %>%
  as_tibble() %>% 
  group_by(Sample, OS, OStime) %>%
  dplyr::select(Sample,OS,OStime) %>%
  distinct() %>%
  mutate(lacksSE = Sample %in% e4less_samples)
survival_data$OSSurv = Surv(time=survival_data$OStime, event=survival_data$OS)
fit <- survfit(OSSurv ~ lacksSE, data = survival_data)
ggsurvplot(fit, data = survival_data, pval = TRUE, palette="Set1", risk.table = T, legend.title = "", ylab="Overall Survival", xlab = "Days")


fit.coxph <- coxph(OSSurv ~ ALKCoAmpl, data = survival_data)
ggforest(fit.coxph, data = survival_data)

mna_profiles_allchr[survival_data %>% filter(nFragment > 3) %>% .$Sample] %>% unlist() %>% as_tibble %>% View


# ------------------------------------------------------------------------------
# Fold Change over background
# ------------------------------------------------------------------------------

fc_bw = read_bigwig("/Volumes/Elements/MYCNAmplicon/Results/CopyNumber_FCOverChance.bigwig") 

plot_xmin = 0
plot_xmax = 40000000


amplicons.fig = 
  mna_profiles_allchr_df %>%
  full_join(fragments_by_sample, by="Sample") %>% 
  mutate(containsE4 = Sample %in% e4less_samples) %>%
  mutate(Sample = forcats::fct_reorder(Sample, -nFragments)) %>% 
  filter(seqnames == "chr2") %>% 
  filter(start<50000000) %>% 
  ggplot(aes(x=start, y=Sample)) + 
  geom_errorbarh(aes(xmin=start, xmax=end, ymin=Sample, ymax=Sample), height=0) +
  theme_kons2() + 
  xlim(plot_xmin,plot_xmax) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y=element_blank()) +
  
  fc.fig = fc_bw %>%
  as_tibble() %>%
  mutate(score_pos = pmax(score, 1),
         score_neg = pmin(score, 1)) %>% 
  ggplot(aes(x = start)) + 
  geom_ribbon(aes(ymin=1, ymax=score_pos), fill="firebrick3") +
  geom_ribbon(aes(ymin=score_neg, ymax=1), fill="steelblue") +
  xlim(plot_xmin,plot_xmax) +
  theme_kons2()

annot_regions_df = data.frame(
  "name" = c("ODC1", "MYCN", "ALK"),
  "start" = c(10580094, 16080683, 29415640),
  "end" = c(10588630, 16087129, 30144477),
  "this_color" = c("black", "black", "black"),
  "this_fill" = c("black", "black", "black")
)
annot_regions_df$dummy = 1

annot_regions.fig = 
  annot_regions_df %>%
  ggplot(aes(x=start, y=dummy)) +
  geom_rect(xmin=annot_regions_df$start, xmax=annot_regions_df$end, ymin=-Inf, ymax=Inf, color=annot_regions_df$this_color, fill=annot_regions_df$this_fill) + 
  theme_kons2() +
  theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank()) +
  xlim(plot_xmin,plot_xmax) 

se.df = 
  read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_lowMYCN_CRCdriven_SE.bed") %>% 
  as_tibble() %>% 
  mutate(dummy=1, this_color="gold", this_fill="gold") %>% 
  filter(seqnames == "chr2")

se.fig = se.df %>%
  ggplot(aes(x=start, y=dummy)) +
  geom_rect(xmin=se.df$start, xmax=se.df$end, ymin=-Inf, ymax=Inf, color=NA, fill=se.df$this_fill) + 
  theme_kons2() +
  theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank()) +
  xlim(plot_xmin,plot_xmax) 

ggsave("/Volumes/Elements/MYCNAmplicon/Results/FCFigure_allChr2.pdf",
       egg::ggarrange(
         annot_regions.fig+
           theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
         se.fig+
           theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
         fc.fig,
         nrow = 3, 
         heights = c(0.2,0.2,1)),
       height=1.5, width=7, onefile = FALSE, useDingbats = F)


ggsave("/Volumes/Elements/MYCNAmplicon/Results/FCFigure_allChr2.pdf",
       egg::ggarrange(
         fc.fig,
         annot_regions.fig+
           theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
         se.fig+
           theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
         nrow = 3, 
         heights = c(1,1,0.2,0.2)),
       height=1.5, width=7, onefile = FALSE, useDingbats = F)

# ------------------------------------------------------------------------------  
# plot all amplicons
# ------------------------------------------------------------------------------

mna_profiles_allchr_df = mna_profiles_allchr %>% unlist()
mna_profiles_allchr_df = unname(mna_profiles_allchr_df)
mna_profiles_allchr_df = as_tibble(mna_profiles_allchr_df)

fragments_by_sample = mna_profiles_allchr_df %>%
  mutate(width = end-start) %>% 
  group_by(Sample) %>% 
  summarise(nFragments = dplyr::n(),
            FragmentLength = sum(width)) %>% 
  ungroup()

mna_profiles_allchr_df %>%
  full_join(fragments_by_sample, by="Sample") %>% 
  mutate(containsE4 = Sample %in% e4less_samples) %>%
  mutate(Sample = forcats::fct_reorder(Sample, -nFragments)) %>% 
  filter(seqnames == "chr2") %>% 
  filter(start<50000000) %>% 
  ggplot(aes(x=start, y=Sample)) + 
  geom_errorbarh(aes(xmin=start, xmax=end), color="firebrick3", height=0) +
  theme_kons2() + 
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ylab("") + xlab("") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y=element_blank()) +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/Amplicons_Chr2_1_50M.pdf",
         height = 1.5, width = 3)

mna_profiles_allchr_df %>%
  full_join(fragments_by_sample, by="Sample") %>% 
  mutate(containsE4 = Sample %in% e4less_samples) %>%
  mutate(Sample = forcats::fct_reorder(Sample, -nFragments)) %>% 
  filter(seqnames == "chr2") %>% 
  filter(start<50000000) %>% 
  ggplot(aes(x=start, y=Sample)) + 
  geom_errorbarh(aes(xmin=start, xmax=end, ymin=Sample, ymax=Sample), height=0) +
  theme_kons2() + 
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ylab("") + xlab("") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y=element_blank()) +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/Amplicons_Chr2_1_50M_Black.pdf",
         height = 1, width = 3)


mna_profiles_allchr_df %>%
  full_join(fragments_by_sample, by="Sample") %>% 
  mutate(containsE4 = Sample %in% e4less_samples) %>%
  mutate(Sample = forcats::fct_reorder(Sample, -nFragments)) %>% 
  filter(seqnames == "chr2") %>% 
  #filter(start<50000000) %>% 
  ggplot(aes(x=start, y=Sample)) + 
  geom_errorbarh(aes(xmin=start, xmax=end)) +
  theme_kons2() + 
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y=element_blank()) +
  facet_zoom(xlim=c(13000000,20000000)) 


mna_profiles_allchr_df %>%
  full_join(fragments_by_sample, by="Sample") %>% 
  mutate(e4less = Sample %in% e4less_samples) %>%
  mutate(Sample = forcats::fct_reorder(Sample, -e4less)) %>% 
  filter(seqnames == "chr2") %>% 
  filter(start<40000000) %>% 
  ggplot(aes(x=start, y=Sample)) + 
  geom_errorbarh(aes(xmin=start, xmax=end, color=e4less), height=0, size=0.2) +
  theme_kons2() + 
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ylab("") + xlab("") +
  scale_color_manual(values=c("FALSE"="steelblue", "TRUE"="firebrick3")) +
  guides(color=F) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y=element_blank()) +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/Amplicons_Chr2_1_40M_e4less.pdf",
         height = 1.5, width = 3)

mna_profiles_allchr_df %>%
  full_join(fragments_by_sample, by="Sample") %>% 
  mutate(e4less = Sample %in% e4less_samples) %>%
  mutate(Sample = forcats::fct_reorder(Sample, -e4less)) %>% 
  filter(seqnames == "chr2") %>% 
  filter(start<40000000) %>% 
  ggplot(aes(x=start, y=Sample)) + 
  geom_errorbarh(aes(xmin=start, xmax=end), height=0, size=0.2) +
  theme_kons2() + 
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ylab("") + xlab("") +
  guides(color=F) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y=element_blank()) +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/Amplicons_Chr2_1_40M.pdf",
         height = 1.5, width = 3)


# ------------------------------------------------------------------------------
# Generate Gene profiles
# ------------------------------------------------------------------------------

library(biomaRt)
library(dplyr)
library(GenomicRanges)

mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl", 
                   mirror = "useast", 
                   GRCh = 37)
attributes <- c("hgnc_symbol", "start_position", "end_position", "gene_biotype")
filters <- c("chromosome_name","start","end")
which_chr = "2"
interval_start = 1
interval_end = 50000000
genes = getBM(attributes=attributes, filters=filters, values=list(chromosome=as.character(which_chr),start=as.character(interval_start),end=as.character(interval_end)), mart=mart) %>% 
  as_tibble() %>% 
  filter(hgnc_symbol != "", gene_biotype == "protein_coding") %>%
  dplyr::select(hgnc_symbol, start_position, end_position)
colnames(genes) = c("gene", "start", "end")
genes$seqnames = "chr2"
genes= makeGRangesFromDataFrame(genes, keep.extra.columns = T)
genes[genes$gene == "MYCN"]

# ------------------------------------------------------------------------------
# only e4less
# ------------------------------------------------------------------------------

genes_df = 
  data.frame(
    "start" = c(15307032, 15731302, 16060521, 16080683, 16190549, 16730727),
    "end" = c(15701454, 15771235, 16076139, 16087129, 16225923, 16847599),
    "name" = c("NBAS", "DDX1", "MYCNUT", "MYCN", "GACAT3", "FAM49A"),
    "class" = rep("gene", 6)
  )
genes.fig = genes_df %>%
  ggplot(aes(x=start, y=1)) +
  geom_rect(xmin=genes_df$start, xmax=genes_df$end, ymin=-Inf, ymax=Inf, color=NA, fill="black") + 
  theme_kons2() + 
  theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank()) +
  xlim(gene_of_interest_start-500000,gene_of_interest_end + 500000)
cRE = read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_nMNA_MYCN_Enhancers.bed") %>% as_tibble()
cRE.fig = cRE %>%
  ggplot(aes(x=start, y=score)) +
  geom_rect(aes(xmin=start, xmax=end), ymin=-Inf, ymax=Inf, color="firebrick3", fill="firebrick3") + 
  theme_kons2() + 
  theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank()) +
  xlim(gene_of_interest_start-500000,gene_of_interest_end + 500000)
cRE.fig 
cRE$class = "cRE"
cRE_genes_plot_df = 
  bind_rows(cRE, genes_df) %>% 
  mutate(dummy=1) %>% 
  dplyr::select(start, end, class, dummy) %>%
  mutate(this_color = ifelse(class == "gene", NA, "firebrick3"),
         this_fill = ifelse(class == "gene", "black", "firebrick3"))
cRE_genes.fig = 
  cRE_genes_plot_df %>%
  ggplot(aes(x=start, y=dummy)) +
  geom_rect(xmin=cRE_genes_plot_df$start, xmax=cRE_genes_plot_df$end, ymin=-Inf, ymax=Inf, color=cRE_genes_plot_df$this_color, fill=cRE_genes_plot_df$this_fill) + 
  theme_kons2() +
  theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank()) +
  xlim(gene_of_interest_start-500000,gene_of_interest_end + 500000)

real_profiles.fig =
  mna_profiles_binned_tb %>%
  mutate(ise4less = Name %in% e4less_samples) %>%
  group_by(start, ise4less) %>%
  summarise(PercentSamples = 100*mean(isAmp),
            nSamples = sum(isAmp)) %>%
  ungroup() %>%
  ggplot(aes(x=start, y=PercentSamples, color=ise4less)) +
  geom_line() + 
  ylab("Samples [%]") + 
  #xlim(gene_of_interest_start-2000000,gene_of_interest_end + 2000000) +
  theme_kons2() 
real_profiles.fig

ggsave("/Volumes/Elements/MYCNAmplicon/Results/e4Less_Aggregate_Copy_Number_Profile.pdf",
       egg::ggarrange(
         genes.fig+
           theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
         cRE.fig+
           theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
         real_profiles.fig,
         nrow = 3, 
         heights = c(0.05,0.05,1)),
       height=2, width=3, onefile = FALSE, useDingbats = F)


set.seed(42)
e4less_samples_mna_profiles = mna_profiles[e4less_samples]
roi = GRanges(seqnames = "chr2",
              ranges = IRanges(start=16080683, end=16087129))
roi_length = width(roi)
randomized_profiles = list()
reps = 100
nrandomsamples = reps*length(e4less_samples_mna_profiles)
library(parallel)
cores = detectCores()
e4less_randomized_profiles = mclapply(1:nrandomsamples,
                               function (i) {
                                 set.seed(i)
                                 random_patient = e4less_samples_mna_profiles[[sample(names(e4less_samples_mna_profiles), 1)]]
                                 random_patient = random_patient[width(random_patient)>roi_length]
                                 random_interval = random_patient[sample(length(random_patient), 1, prob=width(random_patient)-roi_length)]
                                 random_offset = sample.int(width(random_interval)-roi_length,1)
                                 absolute_offset = start(roi) - (start(random_interval) + random_offset)
                                 random_patient = trim(shift(random_patient, absolute_offset))
                                 return(random_patient)
                               },
                               mc.cores = cores)
e4less_randomized_profiles_cov = lapply(e4less_randomized_profiles, coverage)
e4less_randomized_profiles_binned = lapply(e4less_randomized_profiles_cov, function (this_profile_cov) binnedAverage(bins_focus, numvar = this_profile_cov, varname = "isAmp"))
e4less_randomized_profiles_binned_tb = as_tibble(do.call(c, e4less_randomized_profiles_binned))

# ------------------------------------------------------------------------------
# How many % of distal co-amplicons cover SE?
# ------------------------------------------------------------------------------

library(dplyr)
library(plyranges)
rm(list=ls())
load("/Volumes/Elements/MYCNAmplicon/Results/MNAProfiles.Rdata")
lowmycn_crc_se = read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_lowMYCN_CRCdriven_SE.bed")
lowmycn_crc_se = lowmycn_crc_se %>%
  plyranges::filter_by_non_overlaps(GRanges("chr2:16000000-16500000:*"))

lowmycn_crc_enh = read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_lowMYCNExpr_H3K27ac_AggregatePeaks_CRCFactorPositiveOnly.bed")

mycn = GRanges("chr2:16080683-16087129:*")
lowmycn_crc_se_overlaps = list()
lowmycn_crc_enh_overlaps = list()
for (patient in names(mna_profiles)){
  lowmycn_crc_se_overlaps[[patient]] = 
    mna_profiles[[patient]] %>% 
    #filter_by_non_overlaps(mycn) %>% 
    filter_by_overlaps(lowmycn_crc_se) %>% 
    length()
  
  lowmycn_crc_enh_overlaps[[patient]] = 
    mna_profiles[[patient]] %>% 
    #filter_by_non_overlaps(mycn) %>% 
    filter_by_overlaps(lowmycn_crc_enh) %>% 
    length()
}
lowmycn_crc_se_overlaps = unlist(lowmycn_crc_se_overlaps)
length(lowmycn_crc_se_overlaps)
sum(lowmycn_crc_se_overlaps>0)
mean(lowmycn_crc_se_overlaps>0)

lowmycn_crc_enh_overlaps = unlist(lowmycn_crc_enh_overlaps)
length(lowmycn_crc_enh_overlaps)
sum(lowmycn_crc_enh_overlaps>0)
mean(lowmycn_crc_enh_overlaps>0)

# ------------------------------------------------------------------------------
# 
# ------------------------------------------------------------------------------
library(dplyr)
library(plyranges)
rm(list=ls())
load("/Volumes/Elements/MYCNAmplicon/Results/e4less_MNAProfiles.Rdata")
lowmycn_crc_se = read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_lowMYCN_CRCdriven_SE.bed")
lowmycn_crc_enh = read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_lowMYCNExpr_H3K27ac_AggregatePeaks_CRCFactorPositiveOnly.bed")

lowmycn_crc_se = lowmycn_crc_se %>%
  plyranges::filter_by_non_overlaps(GRanges("chr2:16000000-16500000:*"))

lowmycn_crc_se_overlaps = list()
lowmycn_crc_enh_overlaps = list()
for (patient in names(e4less_mna_profiles)){
  lowmycn_crc_se_overlaps[[patient]] = 
    e4less_mna_profiles[[patient]] %>% 
    filter_by_overlaps(lowmycn_crc_se) %>% 
    length()
  
  lowmycn_crc_enh_overlaps[[patient]] = 
    e4less_mna_profiles[[patient]] %>% 
    filter_by_overlaps(lowmycn_crc_enh) %>% 
    length()
}
lowmycn_crc_enh_overlaps = unlist(lowmycn_crc_enh_overlaps)
lowmycn_crc_se_overlaps = unlist(lowmycn_crc_enh_overlaps)
mean(lowmycn_crc_enh_overlaps>0)
mean(lowmycn_crc_se_overlaps>0)

# ------------------------------------------------------------------------------
# 



mna_profiles_binned_tb %>%
  group_by(start) %>%
  summarise(PercentSamples = 100*mean(isAmp),
            nSamples = sum(isAmp)) %>%
  ungroup() %>%
  filter(start>12500000, start<20000000) %>% 
  ggplot(aes(x=start, y=PercentSamples)) +
  geom_line() + 
  theme_kons1() + 
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ylab("Samples [%]") + 
  xlab("chr2") +
  annotate("rect", xmin = 14933889, xmax=15682339, ymin=-Inf, ymax=Inf, fill="orange", alpha=0.4) +
  annotate("text", x = (14933889+15682339)/2, y=Inf, label="FRA2Ctel", vjust=2, size=2) +
  annotate("rect", xmin = 18523760, xmax=19272388, ymin=-Inf, ymax=Inf, fill="orange", alpha=0.4) +
  annotate("text", x = (18523760+19272388)/2, y=Inf, label="FRA2Ccen", vjust=2, size=2) +
  annotate("rect", xmin =16881267,xmax=19399761, ymin=-Inf, ymax=Inf, fill="red", alpha=0.4) +
  annotate("text", x = (16881267+19399761)/2, y=Inf, label="FRA2C", vjust=2, size=2) + 
ggsave("/Volumes/Elements/MYCNAmplicon/Results/ProfileVsFRA.pdf", 
       height=2, width=3, useDingbats=F)
  
  
