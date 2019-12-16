rm(list=ls())
library(dplyr)
library(plyranges)
library(regioneR) # conda install -c bioconda bioconductor-regioner
library(parallel)
library(ggplot2) 
library(wPerm) # conda install -c r r-wperm

blacklist_fname = "/Volumes/Elements/MYCNAmplicon/Data/hg19-blacklist.v2.bed"
hmcan_path = "/Volumes/Elements/nb-cl-chipseq-results/hmcan/"
MACS2_path = "/Volumes/Elements/nb-cl-chipseq-results/MACS2/"
bam_path = "/Volumes/Elements/nb-cl-chipseq-results/bam/"

results_path = "/Volumes/Elements/MYCNAmplicon/Results/analysePeaks/"

encode_mask = data.table::fread(blacklist_fname) %>%
  as_tibble() %>%
  mutate(seqnames = V1, start = V2, end=V3) %>%
  makeGRangesFromDataFrame()

metadata_chipseq = data.frame(
  CellType = c("GICAN", "SH-EP", "SK-N-AS", "GIMEN", "SK-N-SH", "NB69", "SJNB12",
               "SH-SY5Y", "SJNB1", "SK-N-FI", "CLB-GA", "NB-EBc1",
               "LAN1", "CLB-PE", "SK-N-DZ", "CLB-CAR", "CLB-MA",
               "IMR32", "CHP212", "SJNB8", "TR14", #"SK-N-BE2-C",
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

metadata_chipseq$SampleName = paste0("Boeva_", metadata_chipseq$CellType, "_H3K27ac")

metadata_chipseq_MNA = metadata_chipseq %>% filter(Class == "MNA")

metadata_chipseq_MNA$AmpliconSize = NA

metadata_chipseq_MNA$FRIPAmplicon = NA
metadata_chipseq_MNA$FRIPRandomMean = NA
metadata_chipseq_MNA$FRIPRandomSD = NA
metadata_chipseq_MNA$FRIPEmpiricalP = NA
metadata_chipseq_MNA$NumberOfValidSimulations = NA

metadata_chipseq_MNA$NPeaksAmplicon = NA
metadata_chipseq_MNA$NPeaksRandomMean = NA
metadata_chipseq_MNA$NPeaksRandomSD = NA
metadata_chipseq_MNA$NPeaksEmpiricalP = NA

metadata_chipseq_MNA$MeanPeakFCAmplicon = NA
metadata_chipseq_MNA$MeanPeakFCRandomMean = NA
metadata_chipseq_MNA$MeanPeakFCRandomSD = NA
metadata_chipseq_MNA$MeanPeakFCEmpiricalP = NA

metadata_chipseq_MNA$BackgroundChangeInput = NA
metadata_chipseq_MNA$BackgroundChangeH3K27ac = NA

# Fraction of reads in peaks

simulation_overview = list()
for (sample_idx in 11:nrow(metadata_chipseq_MNA)){
  
  sample_name = metadata_chipseq_MNA[[sample_idx,"SampleName"]]
  
  print(sample_name)
  
  amplified_regions = 
    read.table(
      paste0(hmcan_path, 
             sample_name, "/",
             sample_name, "_hmcan_CNV_profile.txt"), 
      header = T) %>%
    filter(Ratio > 10) 
  
  if (nrow(amplified_regions)==0) next
  
  amplified_regions = amplified_regions %>%
    mutate(End = Start + 99999) %>%
    dplyr::select(Chromosome, Start, End) %>%
    makeGRangesFromDataFrame() %>%
    reduce()
  seqlevelsStyle(amplified_regions) = "UCSC"
  
  amplified_regions_size = sum(width(amplified_regions))
  
  wholegenome_bigwig = 
    read_bigwig(paste0(bam_path, 
                       sample_name, ".trimmed.bwa_hg19.rmdup.filtered.bw"))
  wholegenome_bigwig_background = median(wholegenome_bigwig$score, na.rm=T)
  
  wholegenome_inputbigwig = amplified_regions_bigwig = 
    read_bigwig(paste0(bam_path, 
                       gsub("H3K27ac", "Input", sample_name), ".trimmed.bwa_hg19.rmdup.filtered.bw"))
  wholegenome_inputbigwig_background = median(wholegenome_inputbigwig$score, na.rm=T)
  
  amplified_regions_bigwig = 
    read_bigwig(paste0(bam_path, 
                       sample_name, ".trimmed.bwa_hg19.rmdup.filtered.bw"),
                overlap_ranges = amplified_regions)
  amplified_regions_bigwig_background = median(amplified_regions_bigwig$score)
  
  amplified_regions_inputbigwig = 
    read_bigwig(paste0(bam_path, 
                       gsub("H3K27ac", "Input", sample_name), ".trimmed.bwa_hg19.rmdup.filtered.bw"),
                overlap_ranges = amplified_regions)
  amplified_regions_inputbigwig_background = median(amplified_regions_inputbigwig$score)
  
  amplified_regions_peaks = 
    read_narrowpeaks(paste0(MACS2_path,
                            sample_name, "/",
                            sample_name, "_MACS2_peaks.narrowPeak"),
                     overlap_ranges = amplified_regions)
  
  amplified_regions_bigwig_inpeaks = 
    amplified_regions_bigwig %>%
    filter_by_overlaps(amplified_regions_peaks)
  
  real_frip = sum(amplified_regions_bigwig_inpeaks$score) / 
    sum(amplified_regions_bigwig$score)
  
  real_npeaks = length(amplified_regions_peaks)
  
  amplified_regions_peaks$maxBigWigInPeak = amplified_regions_peaks %>%
    group_by_overlaps(amplified_regions_bigwig_inpeaks) %>%
    summarise(maxBigWigInPeak = max(score.query)) %>%
    .$maxBigWigInPeak
  amplicon_background_signal = median(amplified_regions_bigwig$score, na.rm=T)
  amplified_regions_peaks$maxFCInPeakRelativeToRegionBackground = (amplified_regions_peaks$maxBigWigInPeak+1) / (amplicon_background_signal+1) 
  
  real_meanpeakFC = mean(amplified_regions_peaks$maxFCInPeakRelativeToRegionBackground, na.rm=T)
  real_maxpeakFC = max(amplified_regions_peaks$maxFCInPeakRelativeToRegionBackground, na.rm=T)
  
  this_sample_simulation_overview = list()
  this_sample_simulation_overview[[1]] = data.frame(
    "Sample" = sample_name,
    "RandomOrReal" = "Real",
    "SimulationNumber" = NA,
    "RegionBackground" = amplicon_background_signal,
    "RegionFRIP" = real_frip,
    "maxFCInPeakRelativeToRegionBackground" = amplified_regions_peaks$maxFCInPeakRelativeToRegionBackground
  )
  
  randomized_frip = list()
  randomized_npeaks = list()
  randomized_meanpeakFC = list()
  for (i in 1:30){
    
    set.seed(i)
    
    randomized_amplified_regions = randomizeRegions(amplified_regions, 
                                                    genome="hg19",
                                                    mask = encode_mask)
    randomized_amplified_regions_bigwig = 
      read_bigwig(paste0(bam_path, 
                         sample_name, ".trimmed.bwa_hg19.rmdup.filtered.bw"),
                  overlap_ranges = randomized_amplified_regions)
    randomized_amplified_regions_bigwig_background = median(randomized_amplified_regions_bigwig$score)
    
    randomized_amplified_regions_peaks = 
      read_narrowpeaks(paste0(MACS2_path,
                              sample_name, "/",
                              sample_name, "_MACS2_peaks.narrowPeak"),
                       overlap_ranges = randomized_amplified_regions)
    randomized_amplified_regions_bigwig_inpeaks = 
      randomized_amplified_regions_bigwig %>%
      filter_by_overlaps(randomized_amplified_regions_peaks)
    
    randomized_amplified_regions_peaks$maxBigWigInPeak = randomized_amplified_regions_peaks %>%
      group_by_overlaps(randomized_amplified_regions_bigwig_inpeaks) %>%
      summarise(maxBigWigInPeak = max(score.query)) %>%
      .$maxBigWigInPeak
    randomized_amplicon_background_signal = median(randomized_amplified_regions_bigwig$score, na.rm=T)
    randomized_amplified_regions_peaks$maxFCInPeakRelativeToRegionBackground = (randomized_amplified_regions_peaks$maxBigWigInPeak+1) / (randomized_amplicon_background_signal+1) 
    
    randomized_npeaks[[i]] = length(randomized_amplified_regions_peaks)
    
    randomized_frip[[i]] = sum(randomized_amplified_regions_bigwig_inpeaks$score) / 
      sum(randomized_amplified_regions_bigwig$score)
    
    randomized_meanpeakFC[[i]] =  mean(randomized_amplified_regions_peaks$maxFCInPeakRelativeToRegionBackground, na.rm=T)
    
    nrandompeaks = length(randomized_amplified_regions_peaks)
    this_sample_simulation_overview[[i+1]] = data.frame(
      "Sample" = rep(sample_name, nrandompeaks),
      "RandomOrReal" = rep("Random", nrandompeaks),
      "SimulationNumber" = rep(i, nrandompeaks),
      "RegionBackground" = rep(randomized_amplicon_background_signal,nrandompeaks),
      "RegionFRIP" = rep(randomized_frip[[i]], nrandompeaks),
      "maxFCInPeakRelativeToRegionBackground" = 
        randomized_amplified_regions_peaks$maxFCInPeakRelativeToRegionBackground
    )
    
  }
  randomized_frip = unlist(randomized_frip)
  randomized_npeaks = unlist(randomized_npeaks)
  randomized_meanpeakFC = unlist(randomized_meanpeakFC)
  
  metadata_chipseq_MNA[[sample_idx, "AmpliconSize"]] = amplified_regions_size
  
  metadata_chipseq_MNA[[sample_idx, "FRIPAmplicon"]] = real_frip
  metadata_chipseq_MNA[[sample_idx, "FRIPRandomMean"]] = mean(randomized_frip, na.rm=T)
  metadata_chipseq_MNA[[sample_idx, "FRIPRandomSD"]]= sd(randomized_frip, na.rm=T)
  metadata_chipseq_MNA[[sample_idx, "FRIPEmpiricalP"]]= (min(sum(real_frip>randomized_frip, na.rm=T), sum(real_frip<randomized_frip, na.rm=T)) +1) / (sum(!is.na(randomized_frip))+1) 
  
  metadata_chipseq_MNA[[sample_idx, "NPeaksAmplicon"]] = real_npeaks
  metadata_chipseq_MNA[[sample_idx, "NPeaksRandomMean"]] = mean(randomized_npeaks, na.rm=T)
  metadata_chipseq_MNA[[sample_idx, "NPeaksRandomSD"]]= sd(randomized_npeaks, na.rm=T)
  metadata_chipseq_MNA[[sample_idx, "NPeaksEmpiricalP"]]= (min(sum(real_npeaks>randomized_npeaks, na.rm=T), sum(real_npeaks<randomized_npeaks, na.rm=T)) +1) / (sum(!is.na(randomized_npeaks))+1) 
  
  metadata_chipseq_MNA[[sample_idx, "MeanPeakFCAmplicon"]] = real_meanpeakFC
  metadata_chipseq_MNA[[sample_idx, "MeanPeakFCRandomMean"]] = mean(randomized_meanpeakFC, na.rm=T)
  metadata_chipseq_MNA[[sample_idx, "MeanPeakFCRandomSD"]]= sd(randomized_meanpeakFC, na.rm=T)
  metadata_chipseq_MNA[[sample_idx, "MeanPeakFCEmpiricalP"]]= (min(sum(real_meanpeakFC>randomized_meanpeakFC, na.rm=T), sum(real_meanpeakFC<randomized_meanpeakFC, na.rm=T)) +1) / (sum(!is.na(randomized_meanpeakFC))+1) 
  
  metadata_chipseq_MNA[[sample_idx, "BackgroundChangeInput"]] = amplified_regions_inputbigwig_background / wholegenome_inputbigwig_background
  metadata_chipseq_MNA[[sample_idx, "BackgroundChangeH3K27ac"]] = amplified_regions_bigwig_background / wholegenome_bigwig_background
  
  metadata_chipseq_MNA[[sample_idx, "NumberOfValidSimulations"]]= sum(!is.na(randomized_frip))
  print(metadata_chipseq_MNA[sample_idx,])
  
  simulation_overview[[sample_idx]] = do.call(rbind, this_sample_simulation_overview)
  
}
simulation_overview = do.call(rbind, simulation_overview)
metadata_chipseq_MNA$PeakShrink = metadata_chipseq_MNA$MeanPeakFCAmplicon / metadata_chipseq_MNA$MeanPeakFCRandomMean

rm(wholegenome_bigwig,wholegenome_inputbigwig,
   amplified_regions_peaks, amplified_regions_bigwig,
   amplified_regions_inputbigwig)
save.image(paste0(results_path, "Boeva_H3K27ac_PeakSizeAnalysis.Rdata"))


# ------------------------------------------------------------------------------
# Plotting of simulations
# ------------------------------------------------------------------------------
source("/Volumes/Elements/MYCNAmplicon/Code/CustomThemes.R")

simulation_overview %>%
  group_by(Sample, RandomOrReal, SimulationNumber) %>%
  summarise(frip = mean(RegionFRIP, na.rm=T)) %>%
  ungroup() %>% 
  mutate(RandomOrReal = factor(RandomOrReal, levels=c("Random", "Real"))) %>% 
  ggplot(aes(x=Sample, y=frip, color=RandomOrReal, fill=RandomOrReal)) +
  #geom_point(alpha=0.1, position = position_jitterdodge()) + 
  geom_boxplot(fill=NA, outlier.alpha=0) +
  theme_kons2() +
  #scale_y_log10() +
  scale_color_manual(values=c("Random" = "gray50", "Real" = "firebrick3")) +
  scale_fill_manual(values=c("Random" = "gray50", "Real" = "firebrick3")) +
  guides(color=F, fill=F) +
  ggsave(paste0(results_path, "BoevaPeakAnalysis_FRIPRandomVsReal.pdf"),
         height=2,width=3)

simulation_overview %>%
  group_by(Sample, RandomOrReal, SimulationNumber) %>%
  summarise(frip = mean(RegionFRIP, na.rm=T)) %>%
  ungroup() %>% 
  group_by(Sample, RandomOrReal) %>%
  summarise(nDatapoints=dplyr::n())

simulation_overview %>%
  mutate(RandomOrReal = factor(RandomOrReal, levels=c("Random", "Real"))) %>% 
  ggplot(aes(x=Sample, y=maxFCInPeakRelativeToRegionBackground, color=RandomOrReal, fill=RandomOrReal)) +
  #geom_point(alpha=0.1, position = position_jitterdodge()) + 
  geom_boxplot(fill=NA, outlier.alpha=0) +
  theme_kons2() +
  scale_y_log10() +
  scale_color_manual(values=c("Random" = "gray50", "Real" = "firebrick3")) +
  scale_fill_manual(values=c("Random" = "gray50", "Real" = "firebrick3")) +
  guides(color=F, fill=F)+
  ggsave(paste0(results_path, "BoevaPeakAnalysis_PeakFCRandomVsReal.pdf"),
         height=2,width=3)

simulation_overview %>%
  group_by(Sample, RandomOrReal, SimulationNumber) %>%
  summarise(n=dplyr::n()) %>%
  ungroup() %>% 
  mutate(RandomOrReal = factor(RandomOrReal, levels=c("Random", "Real"))) %>% 
  ggplot(aes(x=Sample, y=n, color=RandomOrReal, fill=RandomOrReal)) +
  #geom_point(alpha=0.1, position = position_jitterdodge()) + 
  geom_boxplot(fill=NA, outlier.alpha=0) +
  theme_kons2() +
  scale_color_manual(values=c("Random" = "gray50", "Real" = "firebrick3")) +
  scale_fill_manual(values=c("Random" = "gray50", "Real" = "firebrick3")) +
  guides(color=F, fill=F)+
  ggsave(paste0(results_path, "BoevaPeakAnalysis_PeakNumberRandomVsReal.pdf"),
         height=2,width=3)

simulation_overview %>%
  group_by(Sample, RandomOrReal, SimulationNumber) %>%
  summarise(n=dplyr::n()) %>%
  ungroup() %>% 
  mutate(RandomOrReal = factor(RandomOrReal, levels=c("Random", "Real"))) %>%
  group_by(Sample, RandomOrReal) %>%
  summarise(nDatapoints=dplyr::n())

# ------------------------------------------------------------------------------
# Some random data exploration
# -----------------------------------------------------------------------------

simulation_overview %>%
  mutate(RandomOrReal = factor(RandomOrReal, levels=c("Random", "Real"))) %>%
  group_by(Sample, RandomOrReal) %>%
  summarise(CV_PeakHeight = sd(maxFCInPeakRelativeToRegionBackground, 
                               na.rm=T) / mean(maxFCInPeakRelativeToRegionBackground, 
                                               na.rm=T))

mean(metadata_chipseq_MNA$BackgroundChangeH3K27ac / metadata_chipseq_MNA$BackgroundChangeInput, na.rm=T)

t.test(metadata_chipseq_MNA$BackgroundChangeH3K27ac, metadata_chipseq_MNA$BackgroundChangeInput, paired=T)
# p = 0.09

set.seed(1)
perm.paired.loc(
  metadata_chipseq_MNA$BackgroundChangeH3K27ac[!is.na(metadata_chipseq_MNA$BackgroundChangeH3K27ac)], 
  metadata_chipseq_MNA$BackgroundChangeInput[!is.na(metadata_chipseq_MNA$BackgroundChangeH3K27ac)],
  parameter = mean,
  alternative = "two.sided")
# p = 0.0394

metadata_chipseq_MNA %>%
  ggplot(aes(x = NPeaksAmplicon, y=MeanPeakFCAmplicon)) +
  geom_point()
cor.test(metadata_chipseq_MNA$NPeaksAmplicon, metadata_chipseq_MNA$MeanPeakFCAmplicon)
# cor 0.138

metadata_chipseq_MNA %>%
  ggplot(aes(x = NPeaksAmplicon/AmpliconSize, y=MeanPeakFCAmplicon)) +
  geom_point() +
  scale_y_log10() + 
  scale_x_log10()
cor.test(metadata_chipseq_MNA$NPeaksAmplicon/metadata_chipseq_MNA$AmpliconSize, metadata_chipseq_MNA$MeanPeakFCAmplicon,
         method="spearman")

metadata_chipseq_MNA %>%
  ggplot(aes(x = NPeaksAmplicon/AmpliconSize, y=PeakShrink)) +
  geom_point()
cor.test(metadata_chipseq_MNA$NPeaksAmplicon/metadata_chipseq_MNA$AmpliconSize, metadata_chipseq_MNA$PeakShrink,
         method="spearman")

metadata_chipseq_MNA %>%
  ggplot(aes(x = NPeaksAmplicon/AmpliconSize, y=FRIPAmplicon)) +
  geom_point()
cor.test(metadata_chipseq_MNA$NPeaksAmplicon/metadata_chipseq_MNA$AmpliconSize, metadata_chipseq_MNA$FRIPAmplicon)
# cor 0.76, p=0.01
