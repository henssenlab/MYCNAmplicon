library(GenomicRanges)
library(dplyr)
library(plyranges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Rsamtools)
library(IRanges)
library(regioneR)
library(ggplot2)
source("/Volumes/Elements/MYCNAmplicon/Code/CustomThemes.R")
library(ggforce)

standardchrs = standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)[1:22]
# download from https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg19-blacklist.v2.bed.gz
blacklist_fname = "~/Downloads/hg19-blacklist.v2.bed"
encode_mask = 
  data.table::fread(blacklist_fname) %>%
  as_tibble() %>%
  mutate(seqnames = V1, start = V2, end=V3) %>%
  makeGRangesFromDataFrame()
encode_whitelist = gaps(encode_mask, start=1L, end=seqlengths(BSgenome.Hsapiens.UCSC.hg19))
encode_whitelist = encode_whitelist[strand(encode_whitelist ) == "*"]

# ------------------------------------------------------------------------------
# Define Samples
# ------------------------------------------------------------------------------

# Define Files 
metadata = list()

metadata[["CHP-212"]] = 
  data.frame(
    "amplicon_fname" = "/Volumes/Elements/nb-cl-wgs/nb-cl-wgs-freec/CHP-212/CHP-212_IlluminaWGS.hg19.rmdup.fixRG.bam_CNVs.p.value.MYCNAmplicon.txt",
    "wgs_bw_fname" = "/Volumes/Elements/nb-cl-wgs/nb-cl-wgs-bam/CHP-212_IlluminaWGS.hg19.rmdup.fixRG.bw",
    "atac_bw_fname" = "/Volumes/Elements/nb-cl-atacseq-results/bam/Mundlos_CHP212_ATAC.trimmed.hg19.rmdup.filterednormed.bw",
    "atac_bam_fname" = "/Volumes/Elements/nb-cl-atacseq-results/bam/bam/Mundlos_CHP212_ATAC.hg19.rmdup.bam",
    "atac_bai_fname" = "/Volumes/Elements/nb-cl-atacseq-results/bam/bam/Mundlos_CHP212_ATAC.hg19.rmdup.bai",
    "atac_narrowpeaks_fname" = "/Volumes/Elements/nb-cl-atacseq-results/MACS2/Mundlos_CHP212_ATAC/Mundlos_CHP212_ATAC_MACS2_peaks.narrowPeak",
    "samplename" = "CHP-212")

metadata[["IMR-5-75"]] = 
  data.frame(
    "amplicon_fname" = "/Volumes/Elements/nb-cl-wgs/nb-cl-wgs-freec/IMR-5-75/IMR-5-75_IlluminaWGS.hg19.rmdup.fixRG.bam_CNVs.p.value.MYCNAmplicon.txt",
    "wgs_bw_fname" = "/Volumes/Elements/nb-cl-wgs/nb-cl-wgs-bam/IMR-5-75_IlluminaWGS.hg19.rmdup.fixRG.bw",
    "atac_bw_fname" = "/Volumes/Elements/nb-cl-atacseq-results/bam/Mundlos_IMR575_ATAC.trimmed.hg19.rmdup.filterednormed.bw",
    "atac_bam_fname" = "/Volumes/Elements/nb-cl-atacseq-results/bam/bam/Mundlos_IMR575_ATAC.hg19.rmdup.bam",
    "atac_bai_fname" = "/Volumes/Elements/nb-cl-atacseq-results/bam/bam/Mundlos_IMR575_ATAC.hg19.rmdup.bai",
    "atac_narrowpeaks_fname" = "/Volumes/Elements/nb-cl-atacseq-results/MACS2/Mundlos_IMR575_ATAC/Mundlos_IMR575_ATAC_MACS2_peaks.narrowPeak",
    "samplename" = "IMR-5-75")

metadata[["Kelly"]] = 
  data.frame(
    "amplicon_fname" = "/Volumes/Elements/nb-cl-wgs/nb-cl-wgs-freec/Kelly/Kelly_IlluminaWGS.hg19.rmdup.fixRG.bam_CNVs.p.value.MYCNAmplicon.txt",
    "wgs_bw_fname" = "/Volumes/Elements/nb-cl-wgs/nb-cl-wgs-bam/Kelly_IlluminaWGS.hg19.rmdup.fixRG.bw",
    "atac_bw_fname" = "/Volumes/Elements/nb-cl-atacseq-results/bam/Mundlos_KELLY_ATAC.trimmed.hg19.rmdup.filterednormed.bw",
    "atac_bam_fname" = "/Volumes/Elements/nb-cl-atacseq-results/bam/bam/Mundlos_KELLY_ATAC.hg19.rmdup.bam",
    "atac_bai_fname" = "/Volumes/Elements/nb-cl-atacseq-results/bam/bam/Mundlos_KELLY_ATAC.hg19.rmdup.bai",
    "atac_narrowpeaks_fname" = "/Volumes/Elements/nb-cl-atacseq-results/MACS2/Mundlos_KELLY_ATAC/Mundlos_KELLY_ATAC_MACS2_peaks.narrowPeak",
    "samplename" = "Kelly")

metadata[["LAN-1"]] = 
  data.frame(
    "amplicon_fname" = "/Volumes/Elements/nb-cl-wgs/nb-cl-wgs-freec/LAN-1/LAN-1_IlluminaWGS.hg19.rmdup.fixRG.bam_CNVs.p.value.MYCNAmplicon.txt",
    "wgs_bw_fname" = "/Volumes/Elements/nb-cl-wgs/nb-cl-wgs-bam/LAN-1_IlluminaWGS.hg19.rmdup.fixRG.bw",
    "atac_bw_fname" = "/Volumes/Elements/nb-cl-atacseq-results/bam/Mundlos_LAN1_ATAC.trimmed.hg19.rmdup.filterednormed.bw",
    "atac_bam_fname" = "/Volumes/Elements/nb-cl-atacseq-results/bam/bam/Mundlos_LAN1_ATAC.hg19.rmdup.bam",
    "atac_bai_fname" = "/Volumes/Elements/nb-cl-atacseq-results/bam/bam/Mundlos_LAN1_ATAC.hg19.rmdup.bai",
    "atac_narrowpeaks_fname" = "/Volumes/Elements/nb-cl-atacseq-results/MACS2/Mundlos_LAN1_ATAC/Mundlos_LAN1_ATAC_MACS2_peaks.narrowPeak",
    "samplename" = "LAN-1")

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
#
# Peak height and number analysis
#
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

atacseq_analyis_summary = data.frame(
  "Sample" = names(metadata)
) %>% as_tibble()
atacseq_analyis_summary$AmpliconSize = NA
atacseq_analyis_summary$FRIPAmplicon = NA
atacseq_analyis_summary$FRIPRandomMean = NA
atacseq_analyis_summary$FRIPRandomSD= NA
atacseq_analyis_summary$FRIPEmpiricalP= NA
atacseq_analyis_summary$NPeaksAmplicon = NA
atacseq_analyis_summary$NPeaksRandomMean = NA
atacseq_analyis_summary$NPeaksRandomSD= NA
atacseq_analyis_summary$NPeaksEmpiricalP= NA
atacseq_analyis_summary$MeanPeakFCAmplicon = NA
atacseq_analyis_summary$MeanPeakFCRandomMean = NA
atacseq_analyis_summary$MeanPeakFCRandomSD= NA
atacseq_analyis_summary$MeanPeakFCEmpiricalP= NA
atacseq_analyis_summary$NumberOfValidSimulations= NA

simulation_overview = list()
for (sample_idx in 1:length(metadata)){
  
  this_sample = names(metadata)[sample_idx]
  sample_name = metadata[[sample_idx]]$samplename
  print(sample_name)
  
  amplicon_fname = metadata[[this_sample]]$amplicon_fname %>% as.character()
  wgs_bw_fname = metadata[[this_sample]]$wgs_bw_fname %>% as.character()
  atac_bw_fname = metadata[[this_sample]]$atac_bw_fname %>% as.character()
  atac_bam_fname = metadata[[this_sample]]$atac_bam_fname %>% as.character()
  atac_bai_fname = metadata[[this_sample]]$atac_bai_fname %>% as.character()
  atac_narrowpeaks_fname = metadata[[this_sample]]$atac_narrowpeaks_fname %>% as.character()
  samplename = metadata[[this_sample]]$samplename %>% as.character()
  
  amplicon = read.table(amplicon_fname, sep="\t", header=T)
  amplicon$chr = paste0("chr", amplicon$chr)
  amplicon = makeGRangesFromDataFrame(amplicon, keep.extra.columns = T, ignore.strand=T)
  seqlevels(amplicon , pruning.mode="coarse") = standardchrs
  amplicon = join_overlap_intersect(amplicon, encode_whitelist)
  seqlevelsStyle(amplicon) = "UCSC"
  amplicon_size = sum(width(amplicon))

  # wholegenome_bigwig = 
  #   read_bigwig(paste0(bam_path, 
  #                      sample_name, ".trimmed.bwa_hg19.rmdup.filtered.bw"))
  # wholegenome_bigwig_background = median(wholegenome_bigwig$score, na.rm=T)
  # 
  # wholegenome_inputbigwig = amplified_regions_bigwig = 
  #   read_bigwig(paste0(bam_path, 
  #                      gsub("H3K27ac", "Input", sample_name), ".trimmed.bwa_hg19.rmdup.filtered.bw"))
  # wholegenome_inputbigwig_background = median(wholegenome_inputbigwig$score, na.rm=T)
  
  bins = 
    tileGenome(seqlengths(BSgenome.Hsapiens.UCSC.hg19),
                    tilewidth=50, cut.last.tile.in.chrom = T) 
  bins_amplicon = bins %>%
    filter_by_overlaps(amplicon)
  seqlevels(bins_amplicon, pruning.mode="coarse") = standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)
  
  amplicon_bigwig = read_bigwig(atac_bw_fname, overlap_ranges = amplicon)
  seqlevels(amplicon_bigwig, pruning.mode="coarse") = standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)

  amplicon_bigwig_binned = binnedAverage(bins_amplicon, 
                numvar = coverage(amplicon_bigwig, weight=amplicon_bigwig$score), 
                varname = "score")

  amplicon_bigwig_background = median(amplicon_bigwig$score)
  amplicon_bigwig_background
  amplicon_bigwig_background = median(amplicon_bigwig_binned$score)
  amplicon_bigwig_background
  
  # amplified_regions_inputbigwig = 
  #   read_bigwig(paste0(bam_path, 
  #                      gsub("H3K27ac", "Input", sample_name), ".trimmed.bwa_hg19.rmdup.filtered.bw"),
  #               overlap_ranges = amplified_regions)
  # amplified_regions_inputbigwig_background = median(amplified_regions_inputbigwig$score)

  amplicon_peaks = read_narrowpeaks(atac_narrowpeaks_fname, overlap_ranges = amplicon)
  amplicon_peaks = amplicon_peaks %>% 
    plyranges::select(-name, -score, -signalValue, -pValue, -qValue, -peak) %>%
    reduce()
  
  amplicon_bigwig_inpeaks = 
    amplicon_bigwig %>%
     filter_by_overlaps(amplicon_peaks)
  
  # Get reads on amplicon
  flag = scanBamFlag(isDuplicate = FALSE, isProperPair = TRUE, isPaired = TRUE,
                     isSecondaryAlignment = FALSE, isSupplementaryAlignment = FALSE)
  param <- ScanBamParam(which = amplicon, flag = flag)
  number_of_amplicon_reads = countBam(atac_bam_fname, param = param, index=atac_bai_fname)
  number_of_amplicon_reads = sum(number_of_amplicon_reads$records)
  
  # Get reads in peaks
  param <- ScanBamParam(which = amplicon_peaks, flag = flag)
  number_of_amplicon_peaks_reads = countBam(atac_bam_fname, param = param, index=atac_bai_fname)
  number_of_amplicon_peaks_reads = sum(number_of_amplicon_peaks_reads$records)
  
  # Calculate FRIP
  real_frip = number_of_amplicon_peaks_reads / number_of_amplicon_reads
  
  # Calculate number of peaks
  real_npeaks = length(amplicon_peaks)
  
  # Calculate max bigwig signal in peaks (absolute and relative to background)
  amplicon_peaks$maxBigWigInPeak = 
    amplicon_peaks %>%
    group_by_overlaps(amplicon_bigwig_inpeaks) %>%
    summarise(maxBigWigInPeak = max(score)) %>%
    .$maxBigWigInPeak
  amplicon_background_signal = median(amplicon_bigwig$score, na.rm=T)
  amplicon_peaks$maxFCInPeakRelativeToRegionBackground = (amplicon_peaks$maxBigWigInPeak+1) / (amplicon_background_signal+1) 
  real_meanpeakFC = mean(amplicon_peaks$maxFCInPeakRelativeToRegionBackground, na.rm=T)
  #real_maxpeakFC = max(amplified_regions_peaks$maxFCInPeakRelativeToRegionBackground, na.rm=T)
  
  this_sample_simulation_overview = list()
  this_sample_simulation_overview[[1]] = data.frame(
    "Sample" = sample_name,
    "RandomOrReal" = "Real",
    "SimulationNumber" = NA,
    "RegionBackground" = amplicon_background_signal,
    "RegionFRIP" = real_frip,
    "maxFCInPeakRelativeToRegionBackground" = amplicon_peaks$maxFCInPeakRelativeToRegionBackground
  )
  
  randomized_frip = list()
  randomized_npeaks = list()
  randomized_meanpeakFC = list()
  for (i in 1:3){
    
    set.seed(i)
    
    randomized_amplicon = randomizeRegions(amplicon, 
                                           genome="hg19",
                                           mask = c(encode_mask, amplicon),
                                           per.chromosome = T)
    seqlevels(randomized_amplicon, pruning.mode = "coarse") = seqlevels(amplicon)
    
    randomized_amplicon_bigwig = read_bigwig(atac_bw_fname, overlap_ranges = randomized_amplicon)
    randomized_amplicon_bigwig_background = median(randomized_amplicon_bigwig$score)

    randomized_amplicon_peaks = read_narrowpeaks(atac_narrowpeaks_fname, overlap_ranges = randomized_amplicon)
    randomized_amplicon_peaks = randomized_amplicon_peaks %>% 
      plyranges::select(-name, -score, -signalValue, -pValue, -qValue, -peak) %>%
      reduce()
    
    randomized_amplicon_bigwig_inpeaks = 
      randomized_amplicon_bigwig %>%
      filter_by_overlaps(randomized_amplicon_peaks)
    
    # Get reads on amplicon
    flag = scanBamFlag(isDuplicate = FALSE, isProperPair = TRUE, isPaired = TRUE,
                       isSecondaryAlignment = FALSE, isSupplementaryAlignment = FALSE)
    param <- ScanBamParam(which = randomized_amplicon, flag = flag)
    randomized_number_of_amplicon_reads = countBam(atac_bam_fname, param = param, index=atac_bai_fname)
    randomized_number_of_amplicon_reads = sum(randomized_number_of_amplicon_reads$records)
    
    # Get reads in peaks
    param <- ScanBamParam(which = randomized_amplicon_peaks, flag = flag)
    randomized_number_of_amplicon_peaks_reads = countBam(atac_bam_fname, param = param, index=atac_bai_fname)
    randomized_number_of_amplicon_peaks_reads = sum(randomized_number_of_amplicon_peaks_reads$records)
    
    # Calculate FRIP
    randomized_frip[[i]] = randomized_number_of_amplicon_peaks_reads / randomized_number_of_amplicon_reads
    
    # Calculate number of peaks
    randomized_npeaks[[i]] = length(randomized_amplicon_peaks)
    
    # Calculate max bigwig signal in peaks (absolute and relative to background)
    randomized_amplicon_peaks$maxBigWigInPeak = 
      randomized_amplicon_peaks %>%
      group_by_overlaps(randomized_amplicon_bigwig_inpeaks) %>%
      summarise(maxBigWigInPeak = max(score)) %>%
      .$maxBigWigInPeak
    randomized_amplicon_background_signal = median(randomized_amplicon_bigwig$score, na.rm=T)
    randomized_amplicon_peaks$maxFCInPeakRelativeToRegionBackground = (randomized_amplicon_peaks$maxBigWigInPeak+1) / (randomized_amplicon_background_signal+1) 

    randomized_meanpeakFC[[i]] = mean(randomized_amplicon_peaks$maxFCInPeakRelativeToRegionBackground, na.rm=T)

    nrandompeaks = length(randomized_amplicon_peaks)
    this_sample_simulation_overview[[i+1]] = data.frame(
      "Sample" = rep(sample_name, nrandompeaks),
      "RandomOrReal" = rep("Random", nrandompeaks),
      "SimulationNumber" = rep(i, nrandompeaks),
      "RegionBackground" = rep(randomized_amplicon_background_signal,nrandompeaks),
      "RegionFRIP" = rep(randomized_frip[[i]], nrandompeaks),
      "maxFCInPeakRelativeToRegionBackground" =
        randomized_amplicon_peaks$maxFCInPeakRelativeToRegionBackground
    )

  }
  randomized_frip = unlist(randomized_frip)
  randomized_npeaks = unlist(randomized_npeaks)
  randomized_meanpeakFC = unlist(randomized_meanpeakFC)
  
  atacseq_analyis_summary[[sample_idx, "AmpliconSize"]] = amplicon_size
  atacseq_analyis_summary[[sample_idx, "FRIPAmplicon"]] = real_frip
  atacseq_analyis_summary[[sample_idx, "FRIPRandomMean"]] = mean(randomized_frip, na.rm=T)
  atacseq_analyis_summary[[sample_idx, "FRIPRandomSD"]]= sd(randomized_frip, na.rm=T)
  atacseq_analyis_summary[[sample_idx, "FRIPEmpiricalP"]]= (min(sum(real_frip>randomized_frip, na.rm=T), sum(real_frip<randomized_frip, na.rm=T)) +1) / (sum(!is.na(randomized_frip))+1)

  atacseq_analyis_summary[[sample_idx, "NPeaksAmplicon"]] = real_npeaks
  atacseq_analyis_summary[[sample_idx, "NPeaksRandomMean"]] = mean(randomized_npeaks, na.rm=T)
  atacseq_analyis_summary[[sample_idx, "NPeaksRandomSD"]]= sd(randomized_npeaks, na.rm=T)
  atacseq_analyis_summary[[sample_idx, "NPeaksEmpiricalP"]]= (min(sum(real_npeaks>randomized_npeaks, na.rm=T), sum(real_npeaks<randomized_npeaks, na.rm=T)) +1) / (sum(!is.na(randomized_npeaks))+1)

  atacseq_analyis_summary[[sample_idx, "MeanPeakFCAmplicon"]] = real_meanpeakFC
  atacseq_analyis_summary[[sample_idx, "MeanPeakFCRandomMean"]] = mean(randomized_meanpeakFC, na.rm=T)
  atacseq_analyis_summary[[sample_idx, "MeanPeakFCRandomSD"]]= sd(randomized_meanpeakFC, na.rm=T)
  atacseq_analyis_summary[[sample_idx, "MeanPeakFCEmpiricalP"]]= (min(sum(real_meanpeakFC>randomized_meanpeakFC, na.rm=T), sum(real_meanpeakFC<randomized_meanpeakFC, na.rm=T)) +1) / (sum(!is.na(randomized_meanpeakFC))+1)

  atacseq_analyis_summary[[sample_idx, "NumberOfValidSimulations"]]= sum(!is.na(randomized_frip))
  print(atacseq_analyis_summary[sample_idx,])

  simulation_overview[[sample_idx]] = do.call(rbind, this_sample_simulation_overview)
  
}
simulation_overview = do.call(rbind, simulation_overview)
atacseq_analyis_summary$PeakShrink = atacseq_analyis_summary$MeanPeakFCAmplicon / atacseq_analyis_summary$MeanPeakFCRandomMean
atacseq_analyis_summary %>% View

results_path = "/Volumes/Elements/MYCNAmplicon/Results/AmpliconATACAnalysis/"

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
  ggsave(paste0(results_path, "ATAC_FRIPRandomVsReal.pdf"),
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
  ggsave(paste0(results_path, "ATAC_PeakFCRandomVsReal.pdf"),
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
  ggsave(paste0(results_path, "ATAC_PeakNumberRandomVsReal.pdf"),
         height=2,width=3)

simulation_overview %>%
  group_by(Sample, RandomOrReal, SimulationNumber) %>%
  summarise(n=dplyr::n()) %>%
  ungroup() %>% 
  mutate(RandomOrReal = factor(RandomOrReal, levels=c("Random", "Real"))) %>%
  group_by(Sample, RandomOrReal) %>%
  summarise(nDatapoints=dplyr::n())



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
#
# Fragment size analysis
#
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

for (this_sample in names(metadata)){
  
  amplicon_fname = metadata[[this_sample]]$amplicon_fname %>% as.character()
  wgs_bw_fname = metadata[[this_sample]]$wgs_bw_fname %>% as.character()
  atac_bw_fname = metadata[[this_sample]]$atac_bw_fname %>% as.character()
  atac_bam_fname = metadata[[this_sample]]$atac_bam_fname %>% as.character()
  atac_bai_fname = metadata[[this_sample]]$atac_bai_fname %>% as.character()
  samplename = metadata[[this_sample]]$samplename %>% as.character()
  
  amplicon = read.table(amplicon_fname, sep="\t", header=T)
  amplicon$chr = paste0("chr", amplicon$chr)
  amplicon = makeGRangesFromDataFrame(amplicon, keep.extra.columns = T, ignore.strand=T)
  seqlevels(amplicon , pruning.mode="coarse") = standardchrs
  
  amplicon = join_overlap_intersect(amplicon, encode_whitelist)
  
  # wgs_bw = read_bigwig(wgs_bw_fname)
  # seqlevels(wgs_bw , pruning.mode="coarse") = standardchrs
  # 
  # atac_bw = read_bigwig(atac_bw_fname)
  # seqlevels(atac_bw , pruning.mode="coarse") = standardchrs
  
  # cutoff
  
  # insert sizes on amplicon
  what <- c("isize")
  flag = scanBamFlag(isFirstMateRead = TRUE,
                     isDuplicate = FALSE, isPaired = TRUE,
                     isSecondaryAlignment = FALSE, isSupplementaryAlignment = FALSE)
  param <- ScanBamParam(which = amplicon, what = what, flag = flag)
  amplicon_reads = scanBam(atac_bam_fname, param = param, index=atac_bai_fname)
  names(amplicon_reads) = paste0("ID", 1:length(amplicon_reads))
  amplicon_isizes = abs(unname(unlist(amplicon_reads)))
  amplicon_isizes = data.frame(
    "ISizes" = amplicon_isizes
  )
  amplicon_isizes$RandomizationIndex = 0
  amplicon_isizes$RandomOrReal = "Real"
  
  amplicon_isizes %>%
    filter(!is.na(ISizes), ISizes < 2500) %>% 
    ggplot(aes(x=ISizes)) + 
    geom_density(bw="SJ")
  
  # insert sizes on anti-amplicon
  anti_amplicon = gaps(amplicon, start=1L, end=seqlengths(BSgenome.Hsapiens.UCSC.hg19))
  anti_amplicon = anti_amplicon[strand(anti_amplicon) == "*"]
  anti_amplicon = anti_amplicon[seqnames(anti_amplicon) == "chr2"]
  anti_amplicon = join_overlap_intersect(anti_amplicon, encode_whitelist)
  
  what <- c("isize")
  flag = scanBamFlag(isFirstMateRead = TRUE,
                     isDuplicate = FALSE, isPaired = TRUE,
                     isSecondaryAlignment = FALSE, isSupplementaryAlignment = FALSE)
  param <- ScanBamParam(which = anti_amplicon, what = what, flag = flag)
  anti_amplicon_reads = scanBam(atac_bam_fname, param = param, index=atac_bai_fname)
  names(anti_amplicon_reads) = paste0("ID", 1:length(anti_amplicon_reads))
  anti_amplicon_isizes = abs(unname(unlist(anti_amplicon_reads)))
  anti_amplicon_isizes = data.frame(
    "ISizes" = anti_amplicon_isizes
  )
  anti_amplicon_isizes$RandomizationIndex = -1
  anti_amplicon_isizes$RandomOrReal = "AntiAmplicon"
  anti_amplicon_isizes %>%
    filter(!is.na(ISizes), ISizes < 2500) %>% 
    ggplot(aes(x=ISizes)) + 
    geom_density()
  
  # randomized region
  set.seed(42)
  n_randomizations = 100
  all_randomized_amplicon_isizes = list()
  for (i in 1:n_randomizations){
    randomized_amplicon = circularRandomizeRegions(
      amplicon, 
      genome="hg19",
      mask = c(encode_mask, amplicon))
    what <- c("isize")
    flag = scanBamFlag(isFirstMateRead = TRUE,
                       isDuplicate = FALSE, isPaired = TRUE,
                       isSecondaryAlignment = FALSE, isSupplementaryAlignment = FALSE)
    param <- ScanBamParam(which = randomized_amplicon, what = what, flag = flag)
    randomized_amplicon_reads = scanBam(atac_bam_fname, param = param, index=atac_bai_fname)
    names(randomized_amplicon_reads) = paste0("ID", 1:length(randomized_amplicon_reads))
    randomized_amplicon_isizes = abs(unname(unlist(randomized_amplicon_reads)))
    randomized_amplicon_isizes = data.frame(
      "ISizes" = randomized_amplicon_isizes
    )
    randomized_amplicon_isizes$RandomizationIndex = i
    randomized_amplicon_isizes$RandomOrReal = "Random"
    
    randomized_amplicon_isizes = randomized_amplicon_isizes %>%
      filter(!is.na(ISizes), ISizes<2500)
    
    all_randomized_amplicon_isizes[[i]] = randomized_amplicon_isizes
  }
  all_randomized_amplicon_isizes =  do.call(rbind, all_randomized_amplicon_isizes)
  
  percent_large_fragment_df = 
    rbind(amplicon_isizes, all_randomized_amplicon_isizes, anti_amplicon_isizes) %>%
    filter(!is.na(ISizes), ISizes<2500) %>% 
    group_by(RandomizationIndex, RandomOrReal) %>% 
    summarise(PercentLarge = 100*mean(ISizes>1200, na.rm=T)) 

  percent_large_fragment_df %>% 
    ggplot(aes(x=PercentLarge)) + 
    geom_vline(xintercept = percent_large_fragment_df %>% filter(RandomOrReal == "Real") %>% .$PercentLarge, color="red")+
    geom_vline(xintercept = percent_large_fragment_df %>% filter(RandomOrReal == "AntiAmplicon") %>% .$PercentLarge, color="blue")+
    geom_density()
  
  rbind(amplicon_isizes, 
        all_randomized_amplicon_isizes, 
        anti_amplicon_isizes) %>%
    filter(!is.na(ISizes), ISizes<2500) %>%
    ggplot(aes(x=ISizes, color=RandomOrReal)) +
    geom_density() +
    theme_kons2() 
  
  rbind(amplicon_isizes, 
        all_randomized_amplicon_isizes) %>%
    filter(!is.na(ISizes), ISizes<2500) %>%
    ggplot(aes(x=ISizes, color=RandomOrReal)) +
    geom_density() +
    theme_kons2() +
    facet_zoom(xlim = c(1000,2500), ylim = c(0,0.00025), horizontal=F) +
    theme(strip.background = element_rect(fill="grey90", color=NA)) +
    ggtitle(paste0(samplename)) +
    xlab("ATAC-seq Fragment Length [bp]") +
    ggsave(
      paste0("/Volumes/Elements/MYCNAmplicon/Results/AmpliconATACAnalysis/", samplename, "_ATACFragmentLength_AmpliconVsRandom.pdf"),
      height = 3, width = 3, useDingbats=F)
  
  rbind(amplicon_isizes, 
        anti_amplicon_isizes) %>%
    filter(!is.na(ISizes), ISizes<2500) %>%
    ggplot(aes(x=ISizes, color=RandomOrReal)) +
    geom_density() +
    theme_kons2() +
    facet_zoom(xlim = c(1000,2500), ylim = c(0,0.00025), horizontal=F) +
    theme(strip.background = element_rect(fill="grey90", color=NA)) +
    ggtitle(paste0(samplename)) +
    xlab("ATAC-seq Fragment Length [bp]") +
    ggsave(
      paste0("/Volumes/Elements/MYCNAmplicon/Results/AmpliconATACAnalysis/", samplename, "_ATACFragmentLength_AmpliconVsRestChr2.pdf"),
      height = 3, width = 3, useDingbats=F)
  
  rbind(amplicon_isizes, 
        anti_amplicon_isizes,
        all_randomized_amplicon_isizes) %>%
    filter(!is.na(ISizes), ISizes<2500) %>%
    ggplot(aes(x=ISizes, color=RandomOrReal)) +
    geom_density() +
    theme_kons2() +
    facet_zoom(xlim = c(1000,2500), ylim = c(0,0.00025), horizontal=F) +
    theme(strip.background = element_rect(fill="grey90", color=NA)) +
    ggtitle(paste0(samplename)) +
    xlab("ATAC-seq Fragment Length [bp]") +
    ggsave(
      paste0("/Volumes/Elements/MYCNAmplicon/Results/AmpliconATACAnalysis/", samplename, "_ATACFragmentLength_AmpliconVsBoth.pdf"),
      height = 3, width = 3, useDingbats=F)
  
  rbind(amplicon_isizes, 
        all_randomized_amplicon_isizes, 
        anti_amplicon_isizes) %>%
    filter(!is.na(ISizes), ISizes<2500) %>%
    filter(ISizes>1000) %>%
    ggplot(aes(x=ISizes, color=RandomOrReal)) +
    geom_density() +
    theme_kons2() 
  
  rbind(amplicon_isizes, 
        anti_amplicon_isizes) %>%
    filter(!is.na(ISizes), ISizes<2500) %>%
    filter(ISizes>1000) %>%
    ggplot(aes(x=ISizes, color=RandomOrReal)) +
    geom_density() +
    theme_kons2() 
  
}

