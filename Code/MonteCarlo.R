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
mycn_granges = GRanges(
  seqnames = c("chr2"),
  ranges = IRanges(
    start = c(16080683),
    end = c(16087129)
  )
)
arrayDataIsMNA = vector()
for (sample_idx in 1:length(mna_profiles)){
  mna_profiles[[sample_idx]]$Sample = names(mna_profiles)[sample_idx]
  if (length(findOverlaps(mycn_granges, mna_profiles[sample_idx])) > 0) arrayDataIsMNA = c(arrayDataIsMNA, sample_idx)
}
mna_profiles = mna_profiles[arrayDataIsMNA]

save(mna_profiles, file="/Volumes/Elements/MYCNAmplicon/Results/MNAProfiles.Rdata")

# save e4-less
cRE = read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_nMNA_MYCN_Enhancers.bed") 
cRE[cRE$name == "e4"]


########

rm(list=ls())
library(BSgenome.Hsapiens.UCSC.hg19)
library(plyranges)
library(dplyr)

for (i in 1:30){
  
  load("/Volumes/Elements/MYCNAmplicon/Results/MNAProfiles.Rdata")
  
  randomization_index = i
  fname = paste0("/Volumes/Elements/MYCNAmplicon/Results/SimulationResults/", 
                 "Randomization", as.character(randomization_index), ".Rdata")
  
  set.seed(randomization_index)
  
  mycn = GRanges(seqnames = "chr2",
                 ranges = IRanges(start=16080683, end=16087129))
  mycn_length = width(mycn)
  
  randomized_intervals = list()
  for (this_patient in names(mna_profiles)){
    out_of_bounds = T
    while (out_of_bounds){
      this_patient_intervals_original = mna_profiles[[this_patient]]
      this_patient_intervals = mna_profiles[[this_patient]]
      this_patient_intervals = this_patient_intervals[width(this_patient_intervals)>mycn_length]
      random_interval = this_patient_intervals[sample(length(this_patient_intervals), 1, prob=width(this_patient_intervals)-mycn_length)]
      random_offset = sample.int(width(random_interval)-mycn_length,1)
      absolute_offset = start(mycn) - (start(random_interval) + random_offset)
      this_patient_intervals = shift(this_patient_intervals_original, absolute_offset)
      idx <- GenomicRanges:::get_out_of_bound_index(this_patient_intervals)
      if (sum(width(this_patient_intervals)) == sum(width(trim(this_patient_intervals)))) out_of_bounds = F
      this_patient_intervals = trim(this_patient_intervals)
      randomized_intervals[[this_patient]] = this_patient_intervals
    }
  }
  randomized_intervals = GRangesList(randomized_intervals)
  save(randomized_intervals, file=fname)
}


#####

overall_mna_samples = length(mna_profiles) # 240


cRE = read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_nMNA_MYCN_Enhancers.bed") 

cRE = read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_lowMYCN_CRCdriven_SE.bed") 
cRE$name = paste0("SE", 1:length(cRE))

cRE = GRanges(
  seqnames = c(              "chr2",   "chr2"  ),
  ranges = IRanges(start = c(10580094, 29415640),
                   end =   c(10588630, 30144477))
)
cRE$name = c("ODC1", "ALK")

overlap_frequencies = list()
for (i in 1:30){

  randomization_index = i
  fname = paste0("/Volumes/Elements/MYCNAmplicon/Results/SimulationResults/", 
                 "Randomization", as.character(randomization_index), ".Rdata")
  load(fname)
  
  # overlaps
  randomized_intervals_gr = unlist(randomized_intervals)
  randomized_intervals_gr = unname(randomized_intervals_gr)
  hits = findOverlaps(randomized_intervals_gr, cRE)
  randomized_intervals_gr = randomized_intervals_gr[queryHits(hits)]
  randomized_intervals_gr$subjectHit = subjectHits(hits)
  
  overlap_frequencies[[randomization_index]] =
    randomized_intervals_gr %>%
    as_tibble() %>%
    dplyr::select(Sample, subjectHit) %>%
    mutate(subjectHitName = cRE[randomized_intervals_gr$subjectHit]$name) %>%
    group_by(subjectHitName) %>%
    summarise(
      nSamples = n(),
      PercentSamples = 100 * nSamples / overall_mna_samples) %>%
    mutate(RandomizationIndex = randomization_index)
}
overlap_frequencies = do.call(rbind, overlap_frequencies)
overlap_frequencies$RandomOrReal = "Random"
print(overlap_frequencies)

# Real overlap
mna_profiles_gr = unlist(mna_profiles)
mna_profiles_gr = unname(mna_profiles_gr)
hits = findOverlaps(mna_profiles_gr, cRE)
mna_profiles_gr = mna_profiles_gr[queryHits(hits)]
mna_profiles_gr$subjectHit = subjectHits(hits)
real_overlap_frequencies =
  mna_profiles_gr %>%
  as_tibble() %>%
  dplyr::select(Sample, subjectHit) %>%
  mutate(subjectHitName = cRE[mna_profiles_gr$subjectHit]$name) %>%
  group_by(subjectHitName) %>%
  summarise(
    RealnSamples = n(),
    RealPercentSamples = 100 * RealnSamples / overall_mna_samples) %>%
  mutate(RandomizationIndex = NA)
real_overlap_frequencies$RandomOrReal = "Real"
print(real_overlap_frequencies)
#real_overlap_frequencies %>% View

mean_random_overlap_frequencies = 
  overlap_frequencies %>%
  full_join(real_overlap_frequencies %>% dplyr::select(subjectHitName, RealPercentSamples)) %>% 
  group_by(subjectHitName) %>%
  summarise(random_mean = mean(PercentSamples), 
            random_sd = sd(PercentSamples),
            real = RealPercentSamples[[1]],
            fc_real_over_random = (real + 1) / (random_mean + 1), 
            pval = (sum(PercentSamples >= RealPercentSamples) + 1) / (dplyr::n() + 1)
  ) %>% 
  full_join(cRE %>% as_tibble() %>% mutate(subjectHitName = name))
print(mean_random_overlap_frequencies)
View(mean_random_overlap_frequencies)


# ------------------------------------------------------------------------------
# junk
# ------------------------------------------------------------------------------

randomized_intervals_df = randomized_intervals %>% unlist()
randomized_intervals_df = unname(randomized_intervals_df)
randomized_intervals_df = as_tibble(randomized_intervals_df)
randomized_intervals_df %>%
  ggplot(aes(x=start, y=Sample)) + 
  geom_errorbarh(aes(xmin=start, xmax=end), color="firebrick3", height=0) +
  theme_kons2() + 
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ylab("") + xlab("") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y=element_blank()) +
  facet_zoom(xlim=c(16080683, 16087129))

# randomized and mna_profiles have the same numbers of intervals
x = sapply(mna_profiles, length) 
y = sapply(randomized_intervals, length)
sum(x!=y)

# randomized and mna_profiles have the same widthsum of intervals
x = sapply(mna_profiles, function (intervals) sum(width(intervals), na.rm=T)) 
y = sapply(randomized_intervals, function (intervals) sum(width(intervals), na.rm=T))
sum(x!=y)

randomized_intervals_df = randomized_intervals %>% unlist()
randomized_intervals_df = unname(randomized_intervals_df)
randomized_intervals_df = as_tibble(randomized_intervals_df)

# randomized and mna_profiles have the same numbers of intervals
x = sapply(mna_profiles, length) 
y = sapply(randomized_intervals, length)
sum(x!=y)

# randomized and mna_profiles have the same widthsum of intervals
x = sapply(mna_profiles, function (intervals) sum(width(intervals), na.rm=T)) 
y = sapply(randomized_intervals, function (intervals) sum(width(intervals), na.rm=T))
sum(x!=y)

# ------------------------------------------------------------------------------
overlap_frequencies = 
  mclapply(
    1:number_of_randomizations, 
    funtion (i){
      randomization_index = i
      fname = paste0("/fast/users/helmsauk_c/scratch/MYCNAmpliconRandomization/RandomizedProfiles/RandomizedProfiles", as.character(randomization_index), ".Rdata")
      
      load(fname)
      
      # overlaps
      randomized_intervals_gr = unlist(randomized_intervals)
      randomized_intervals_gr = unname(randomized_intervals_gr)
      hits = findOverlaps(randomized_intervals_gr, cRE)
      randomized_intervals_gr = randomized_intervals_gr[queryHits(hits)]
      randomized_intervals_gr$subjectHit = subjectHits(hits)
      
      this_overlap_frequencies =
        randomized_intervals_gr %>%
        as_tibble() %>%
        dplyr::select(Sample, subjectHit) %>%
        mutate(subjectHitName = cRE[randomized_intervals_gr$subjectHit]$name) %>%
        group_by(subjectHitName) %>%
        summarise(
          nSamples = n(),
          PercentSamples = 100 * nSamples / overall_mna_samples) %>%
        mutate(RandomizationIndex = randomization_index)
      
      return(this_overlap_frequencies)
    },
    mc.cores=cores)


# ------------------------------------------------------------------------------

load("/Volumes/Elements/MYCNAmplicon/Results/MNAProfiles.Rdata")

randomization_index = 5
fname = paste0("/Volumes/Elements/MYCNAmplicon/Results/SimulationResults/", 
               "Randomization", as.character(randomization_index), ".Rdata")

set.seed(randomization_index)

# this is a must
mycn = GRanges(seqnames = "chr2",
               ranges = IRanges(start=16080683, end=16087129))
mycn_length = width(mycn)

# this is forbidden
e4 = read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_nMNA_MYCN_Enhancers.bed") %>%
  filter(name == "e4")
e4_length = width(e4)

mycn_e4_start_diff = start(e4) - start(mycn)

this_patient_intervals_original = mna_profiles[[9]]
possible_starts_contain = this_patient_intervals_original
possible_starts_contain = possible_starts_contain[width(possible_starts_contain)>=width(mycn)]
end(possible_starts_contain) = end(possible_starts_contain)-width(mycn)
possible_starts_notcontain = gaps(this_patient_intervals_original)
possible_starts_notcontain = possible_starts_notcontain[strand(possible_starts_notcontain) == "*"] # now I have the gaps
possible_starts_notcontain = possible_starts_notcontain[width(possible_starts_notcontain)>=width(mycn)] # only long enough gaps
end(possible_starts_notcontain) = end(possible_starts_notcontain)-width(e4)
possible_starts_notcontain = trim(shift(possible_starts_notcontain, -mycn_e4_start_diff))
possible_starts = GRanges(coverage(c(possible_starts_contain, possible_starts_notcontain))) %>% filter(score == 2)
random_start_interval = possible_starts[sample(length(possible_starts), 1, prob=width(possible_starts))]
random_mycn_start = start(random_start_interval) + sample.int(width(random_start_interval),1)
absolute_shift = random_mycn_start - start(mycn)
randomized_interval = trim(shift(this_patient_intervals_original, absolute_shift))

set.seed(42)
randomized_intervals = list()
for (i in 1:length(mna_profiles)){
  out_of_bounds = T
  while (out_of_bounds){
    this_patient_intervals_original = mna_profiles[[i]]

    possible_starts_contain = this_patient_intervals_original
    possible_starts_contain = possible_starts_contain[width(possible_starts_contain)>=width(mycn)]
    end(possible_starts_contain) = end(possible_starts_contain)-width(mycn)
    
    possible_starts_notcontain = gaps(this_patient_intervals_original)
    possible_starts_notcontain = possible_starts_notcontain[strand(possible_starts_notcontain) == "*"] # now I have the gaps
    possible_starts_notcontain = possible_starts_notcontain[width(possible_starts_notcontain)>=width(mycn)] # only long enough gaps
    end(possible_starts_notcontain) = end(possible_starts_notcontain)-width(e4)
    possible_starts_notcontain = trim(shift(possible_starts_notcontain, -mycn_e4_start_diff))
    
    possible_starts = GRanges(coverage(c(possible_starts_contain, possible_starts_notcontain))) %>% filter(score == 2)
    random_start_interval = possible_starts[sample(length(possible_starts), 1, prob=width(possible_starts))]
    random_mycn_start = start(random_start_interval) + sample.int(width(random_start_interval),1)
    absolute_shift = start(mycn) - random_mycn_start 
    if (sum(width(this_patient_intervals_original)) == sum(width(trim(shift(this_patient_intervals_original, absolute_shift))))) out_of_bounds = F
    this_patient_intervals = trim(shift(this_patient_intervals_original, absolute_shift))
    randomized_intervals[[i]] = this_patient_intervals
  }
}
randomized_intervals = GRangesList(randomized_intervals)



randomized_intervals_df = randomized_intervals %>% unlist()
randomized_intervals_df = unname(randomized_intervals_df)
#randomized_intervals_df = do.call(c, randomized_intervals_df)
randomized_intervals_df = as_tibble(randomized_intervals_df)
randomized_intervals_df %>%
  ggplot(aes(x=start, y=Sample)) + 
  geom_errorbarh(aes(xmin=start, xmax=end), color="firebrick3", height=0) +
  theme_kons2() + 
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ylab("") + xlab("") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y=element_blank()) +
  #facet_zoom(xlim=c(16080683, 16087129))
  facet_zoom(xlim=c(start(e4), end(e4)))


View(randomized_intervals_df)

# ------------------------------------------------------------------------------

set.seed(42)
mycn = GRanges(seqnames = "chr2",
               ranges = IRanges(start=16080683, end=16087129))
mycn_length = width(mycn)

e4 = GRanges(seqnames = "chr2",
             ranges = IRanges(start=16375277, end=16385344))
e4_length = width(e4)

mycn_e4_start_diff = start(e4) - start(mycn)

randomized_intervals = list()
for (pname in names(mna_profiles)){
  out_of_bounds = T
  while (out_of_bounds){
    this_patient_intervals_original = mna_profiles[[pname]]
    
    possible_starts_contain = this_patient_intervals_original
    possible_starts_contain = possible_starts_contain[width(possible_starts_contain)>=width(mycn)]
    end(possible_starts_contain) = end(possible_starts_contain)-width(mycn)
    
    possible_starts_contain2 = this_patient_intervals_original
    possible_starts_contain2 = possible_starts_contain2[width(possible_starts_contain2)>=width(e4)]
    end(possible_starts_contain2) = end(possible_starts_contain2)-width(e4)
    possible_starts_contain2 = trim(shift(possible_starts_contain2, -mycn_e4_start_diff))
    
    possible_starts = GRanges(coverage(c(possible_starts_contain, possible_starts_contain2))) %>% filter(score == 2)
    
    if (length(possible_starts)==0){
      out_of_bounds = F
      next
    }

    random_start_interval = possible_starts[sample(length(possible_starts), 1, prob=width(possible_starts))]
    random_mycn_start = start(random_start_interval) + sample.int(width(random_start_interval),1)
    absolute_shift = start(mycn) - random_mycn_start
    if (sum(width(this_patient_intervals_original)) == sum(width(trim(shift(this_patient_intervals_original, absolute_shift))))) out_of_bounds = F
    this_patient_intervals = trim(shift(this_patient_intervals_original, absolute_shift))
    randomized_intervals[[pname]] = this_patient_intervals
  }
}
randomized_intervals = GRangesList(randomized_intervals)
randomized_intervals_df = randomized_intervals %>% unlist()
randomized_intervals_df = unname(randomized_intervals_df)
#randomized_intervals_df = do.call(c, randomized_intervals_df)
randomized_intervals_df = as_tibble(randomized_intervals_df)
randomized_intervals_df %>%
  ggplot(aes(x=start, y=Sample)) + 
  geom_errorbarh(aes(xmin=start, xmax=end), color="firebrick3", height=0) +
  theme_kons2() + 
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ylab("") + xlab("") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y=element_blank()) +
  #facet_zoom(xlim=c(16080683, 16087129))
  facet_zoom(xlim=c(start(e4), end(e4)))


