library(plyranges)
library(dplyr)
library(ggplot2)
source("/Volumes/Elements/MYCNAmplicon/Code/CustomThemes.R")
source("/Volumes/Elements/MYCNAmplicon/Code/hasOverlap.R")

blacklist = read.table("/Volumes/Elements/MYCNAmplicon/Data/hg19-blacklist.v2.bed",
                       header = F, sep="\t")
colnames(blacklist) = c("Chr", "Start", "End", "Class")

load("/Volumes/Elements/MYCNAmplicon/Data/UnfilteredSvabaSVs.Rdata")

# ------------------------------------------------------------------------------
# Read fragiles sites
# ------------------------------------------------------------------------------

cfs = read.table("/Volumes/Elements/MYCNAmplicon/Data/chr2_fragilesites_liftOvertohg19.bed",
                 header = F, sep="\t")

# Kumar et al.: FRA2C chr2:16881267-19399761 
# Blumrich et al.: FRA2Ctel chr2:14933889-15682339
# Blumrich et al.: FRAC2cen: chr2:18523760-19272388

# annotate("rect", xmin = 16881267, xmax=19399761, ymin=0, ymax=Inf, fill="blue", alpha=0.2) +
# annotate("text", x = (16881267+19399761)/2, y=Inf, label="FRA2C\nKumar et al.", vjust=2, size=2) +
# annotate("rect", xmin = 14933889, xmax=15682339, ymin=0, ymax=Inf, fill="blue", alpha=0.2) +
# annotate("text", x = (14933889+15682339)/2, y=Inf, label="FRAC2tel\nBlumrich et al.", vjust=2, size=2) +
# annotate("rect", xmin = 18523760, xmax=19272388, ymin=0, ymax=Inf, fill="blue", alpha=0.2) +
# annotate("text", x = (18523760+19272388)/2, y=Inf, label="FRAC2cen\nBlumrich et al.", vjust=2, size=2) +
  
# ------------------------------------------------------------------------------
# MYCN non-amplified NB -- Berlin and Peifer Cohort
# ------------------------------------------------------------------------------

clinical_data = read.table("/Volumes/Elements/MYCNAmplicon/Data/ClinicalData.csv",
                           header=T, sep=",")
mna_samples = unique(as.character(clinical_data[clinical_data$Risk == "MNA", "Sample"]))
non_mna_samples = unique(as.character(clinical_data[clinical_data$Risk != "MNA", "Sample"]))

# source("/Volumes/Elements/MYCNAmplicon/Code/ParseAllSV.Rdata")

svaba_non_mna = svaba %>%
  filter(Sample %in% non_mna_samples) %>% 
  filter(grepl(":1", ID)) #%>% 
  #filter(Filter == "PASS")

breakpoints = 
  data.frame(
    Chr = c(svaba_non_mna$ChrA, svaba_non_mna$ChrB),
    Pos = c(svaba_non_mna$PosA, svaba_non_mna$PosB),
    Sample = c(as.character(svaba_non_mna$Sample), as.character(svaba_non_mna$Sample))
  )

n_samples = length(unique(breakpoints$Sample))
print(n_samples)

breakpoints %>%
  filter(Chr == "chr2", Pos > 15000000, Pos < 17000000) %>%
  mutate(bin = floor(Pos/10000)*10000) %>% 
  group_by(Sample, bin) %>% 
  summarise(nBreakpoints = dplyr::n()) %>% 
  tidyr::complete(Sample, bin, fill=list(nBreakpoints=0)) %>% 
  ggplot(aes(x=bin, y=nBreakpoints, color=Sample)) +
  geom_line() + 
  theme_kons1() + 
  xlab("") + 
  ylab("") +
  guides(color=F)

binsize = 20000
breakpoints %>%
  filter(Chr == "chr2", Pos > 15000000, Pos < 17000000) %>%
  mutate(bin = floor(Pos/binsize)*binsize) %>% 
  group_by(Sample, bin) %>% 
  summarise(nBreakpoints = dplyr::n()) %>% 
  tidyr::complete(Sample, bin, fill=list(nBreakpoints=0)) %>% 
  ggplot(aes(x=bin, y=nBreakpoints)) +
  stat_summary(fun.data=mean_se) + 
  theme_kons1() + 
  xlab("") + 
  ylab("")

breakpoints %>%
  filter(Chr == "chr2", Pos > (gene_of_interest_start-500000), Pos<(gene_of_interest_end + 500000)) %>%
  ggplot(aes(x=Pos, group=Sample)) +
  geom_density() + 
  theme_kons1() + 
  xlab("") + 
  ylab("") + 
  guides(color=F) +
  scale_y_log10()

binsize = 50000
breakpoints %>%
  filter(Chr == "chr2", Pos > 15000000, Pos < 17000000) %>%
  mutate(bin = floor(Pos/binsize)*binsize) %>% 
  group_by(Sample, bin) %>% 
  summarise(nBreakpoints = dplyr::n()) %>% 
  tidyr::complete(Sample, bin, fill=list(nBreakpoints=0)) %>% 
  ungroup() %>% 
  group_by(bin) %>% 
  summarise(meanBreakpoints = mean(nBreakpoints, na.rm=T),
            sdBreakpoints = sd(nBreakpoints, na.rm=T)) %>% 
  ggplot(aes(x=bin, y=meanBreakpoints)) +
  geom_segment(aes(xend=bin+binsize-1, yend=meanBreakpoints)) +
  theme_kons1() + 
  xlab("") + 
  ylab("")

breakpoints %>%
  filter(Chr == "chr2") %>%
  ggplot(aes(x=Pos)) +
  geom_histogram(binwidth = 25000) + 
  theme_kons1() + 
  xlim(15000000,17000000) + 
  xlab("") + 
  ylab("") + 
  # annotate("rect", xmin = 16881267, xmax=19399761, ymin=0, ymax=Inf, fill="blue", alpha=0.2) +
  # annotate("text", x = (16881267+19399761)/2, y=Inf, label="FRA2C\nKumar et al.", vjust=2, size=2) +
  # annotate("rect", xmin = 14933889, xmax=15682339, ymin=0, ymax=Inf, fill="blue", alpha=0.2) +
  # annotate("text", x = (14933889+15682339)/2, y=Inf, label="FRAC2tel\nBlumrich et al.", vjust=2, size=2) +
  # annotate("rect", xmin = 18523760, xmax=19272388, ymin=0, ymax=Inf, fill="blue", alpha=0.2) +
  # annotate("text", x = (18523760+19272388)/2, y=Inf, label="FRAC2cen\nBlumrich et al.", vjust=2, size=2) +
  #ggtitle("Breakpoints around MYCN") + 
  ggsave("/Volumes/Elements/MYCNAmplicon/Figures/NB_nonMNA_BreakpointDistribution.pdf",
         height = 2, width = 5, useDingbats = F)

breakpoints %>%
  filter(Chr == "chr2") %>%
  ggplot(aes(x=Pos)) +
  geom_density(bw=1000) + 
  theme_kons1() + 
  xlim(15000000,17000000) + 
  xlab("") + 
  ylab("") + 
  # annotate("rect", xmin = 16881267, xmax=19399761, ymin=0, ymax=Inf, fill="blue", alpha=0.2) +
  # annotate("text", x = (16881267+19399761)/2, y=Inf, label="FRA2C\nKumar et al.", vjust=2, size=2) +
  # annotate("rect", xmin = 14933889, xmax=15682339, ymin=0, ymax=Inf, fill="blue", alpha=0.2) +
  # annotate("text", x = (14933889+15682339)/2, y=Inf, label="FRAC2tel\nBlumrich et al.", vjust=2, size=2) +
  # annotate("rect", xmin = 18523760, xmax=19272388, ymin=0, ymax=Inf, fill="blue", alpha=0.2) +
  # annotate("text", x = (18523760+19272388)/2, y=Inf, label="FRAC2cen\nBlumrich et al.", vjust=2, size=2) +
  #ggtitle("Breakpoints around MYCN") + 
  ggsave("/Volumes/Elements/MYCNAmplicon/Figures/NB_nonMNA_BreakpointDistribution_Density.pdf",
         height = 2, width = 5, useDingbats = F)


# ------------------------------------------------------------------------------
# Pediatric Pan Cancer Cohort
# ------------------------------------------------------------------------------

# I do not have copy number information that is why I do not know which
# of the tumors are MYCN-amplified

parse_pedpancan_allsv = function(fname){
  d = read.table(fname, header = T, sep=";")
  d$SVCaller = "PedPanCanPipeline"
  d$Cohort = "PedPanCan"
  d$ChrA = as.character(d$ChrA)
  d$ChrB = as.character(d$ChrB)
  d
}
sv = parse_pedpancan_allsv("/Volumes/Elements/MYCNAmplicon/Data/PedPanCanSVs.csv")
meta = read.table("/Volumes/Elements/MYCNAmplicon/Data/PedPanCanMeta.csv", sep=";", header=T)
wgs_samples = meta %>% filter(seq_type == "wgs") %>% .$sample
sv = sv %>% filter(Sample %in% wgs_samples)

breakpoints = 
  data.frame(
    Chr = c(sv$ChrA, sv$ChrB),
    Pos = c(sv$PosA, sv$PosB),
    sample = c(as.character(sv$Sample), as.character(sv$Sample))
) %>%
left_join(meta[, c("sample", "cancer_type", "entity_sub")], by="sample")

n_samples = length(unique(breakpoints$sample))
print(n_samples)
n_cancertypes = length(unique(breakpoints$cancer_type))
print(n_cancertypes)

breakpoints %>%
  filter(Chr == "chr2", entity_sub != "nb_nmyc-amp", cancer_type != "wt", cancer_type != "mb") %>%
  ggplot(aes(x=Pos)) +
  geom_histogram(binwidth = 100000) + 
  theme_kons1() + 
  xlim(10000000,20000000)


  
  
