rm(list=ls())

library(data.table)
library(dplyr)
library(plyranges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)
library(ggforce)
source("/Volumes/Elements/MYCNAmplicon/Code/CustomThemes.R")
source("/Volumes/Elements/MYCNAmplicon/Code/hasOverlap.R")

metadata = data.table::fread("/Volumes/Elements/MYCNAmplicon/Data/mb-arrays.metadata.csv") %>%
  as_tibble() %>%
  mutate(Subgroup = ifelse(Subgroup == "N/A", "Unknown", Subgroup))

# how many medulloblastomas?
data.table::fread("/Volumes/Elements/MYCNAmplicon/Data/mb-arrays.rawcopy.segments.txt") %>%
  as_tibble() %>%
  .$Sample %>%
  unique() %>%
  length()
# 1097 ---> 10 cell lines and 1087 patient samples = samples analysed by MAGIC in Northcott et al.

seg =
  data.table::fread("/Volumes/Elements/MYCNAmplicon/Data/mb-arrays.rawcopy.segments.txt") %>%
  as_tibble() %>%
  dplyr::select(-Allelic.Imbalance) %>%
  filter(Value >= 1.8) %>%
  dplyr::mutate(ID = gsub("^[[:alnum:]_]*-", "", Sample)) %>%
  makeGRangesListFromDataFrame(split.field = "ID", keep.extra.columns=T)
seqlevels(seg, pruning.mode="coarse") = standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)

mycn_granges = GRanges(
  seqnames = c("chr2"),
  ranges = IRanges(
    start = c(16080683),
    end = c(16087129)
  )
)

arrayDataIsMNA = vector()
for (sample_idx in 1:length(seg)){
  if (length(findOverlaps(mycn_granges, seg[sample_idx])) > 0) arrayDataIsMNA = c(arrayDataIsMNA, sample_idx)
}
seg_mna = seg[arrayDataIsMNA]

seg %>%
  unlist() %>% 
  write_bed("/Volumes/Elements/MYCNAmplicon/Results/MB_AllCopyNumberProfiles.bed")
seg_mna %>%
  unlist() %>% 
  write_bed("/Volumes/Elements/MYCNAmplicon/Results/MB_MNANumberProfiles.bed")

seg_mna %>%
  unlist() %>%
  as_tibble() %>% 
  dplyr::select(Sample, seqnames, start, end, Value) %>%
  write.table("/Volumes/Elements/MYCNAmplicon/Results/SourceData/MB_MNA_AmplifiedSegments.txt",
              col.names=T, row.names=F, quote = F, sep = "\t")


bins = tileGenome(seqlengths = seqlengths(BSgenome.Hsapiens.UCSC.hg19),
                  tilewidth = 10000,
                  cut.last.tile.in.chrom = T)
seqlevels(bins, pruning.mode = "coarse") = standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)
seg_mna_cov = lapply(seg_mna, coverage)
seg_mna_binned = lapply(seg_mna_cov, function (this_profile_cov) binnedAverage(bins, numvar = this_profile_cov, varname = "isAmp"))

for (i in 1:length(seg_mna_binned)){
  seg_mna_binned[[i]]$Name = names(seg_mna)[i]
}
seg_mna_binned_tb = do.call(rbind, lapply(seg_mna_binned, as_tibble))
seg_mna_binned_tb$isAmp = seg_mna_binned_tb$isAmp > 0

seg_mna_binned_tb = seg_mna_binned_tb %>%
  inner_join(metadata %>% dplyr::select(ID, Subgroup) %>% dplyr::rename(Name = ID))

seg_mna_binned_tb %>%
  filter(seqnames == "chr2") %>% 
  dplyr::select(Name, Subgroup, seqnames, start, isAmp) %>% 
  mutate(isAmplified = as.logical(isAmp)) %>% 
  dplyr::select(-isAmp) %>% 
  write.table("/Volumes/Elements/MYCNAmplicon/Results/SourceData/MB_MNAsamples_chr2_10kbbins.txt",
            col.names=T, row.names=F, quote = F, sep = "\t")

seg_mna_binned_tb %>%
  .$Name %>%
  unique() %>% 
  length()

seg_mna_binned_tb %>%
  dplyr::select(Name, Subgroup) %>%
  distinct() %>%
  group_by(Subgroup) %>%
  summarise(nMNASamples = n_distinct(Name))

seg_mna_tb = seg_mna %>% unlist() %>% as_tibble()
seg_mna_tb = seg_mna_tb %>% 
  dplyr::mutate(Sample = gsub("^[[:alnum:]_]*-", "", Sample)) %>%
  inner_join(metadata %>% dplyr::select(ID, Subgroup) %>% dplyr::rename(Sample = ID))

# ------------------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------------------
gene_of_interest = "MYCN"
gene_of_interest_ensembl = "ENSG00000134323"
gene_of_interest_chr = "chr2"
gene_of_interest_start = 16080683
gene_of_interest_end = 16087129
window_of_interest = 1250000
view_chr = as.character(seqnames(mycn_granges))
view_start = start(mycn_granges)-window_of_interest
view_end = end(mycn_granges)+window_of_interest
view_granges = GRanges(seqnames = view_chr, ranges = IRanges(start = view_start, end = view_end))

genes_df = 
  data.frame(
    "start" = c(12856998, 14772810, 15307032, 15731302, 16060521, 16080683, 16190549, 16730727, 17691851, 17720393, 17845079, 17935125, 17997763, 18059114, 18735989),
    "end" = c(12882860, 14790933, 15701454, 15771235, 16076139, 16087129, 16225923, 16847599, 17699706, 17838285, 17981509, 17966632, 17998368, 18542882, 18741959),
    "name" = c("TRIB2", "FAM84A", "NBAS", "DDX1", "MYCNUT", "MYCN", "GACAT3", "FAM49A", "RAD51AP2", "VSNL1", "SMC6", "GEN1", "MSGN1", "KCNS3", "RDH14"),
    "class" = rep("gene", 15)
  )

genes.fig = genes_df %>%
  ggplot(aes(x=start, y=1)) +
  geom_rect(xmin=genes_df$start, xmax=genes_df$end, ymin=-Inf, ymax=Inf, color=NA, fill="black") + 
  theme_kons2() + 
  theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank()) +
  xlim(start(view_granges),end(view_granges))

real_profiles.fig =
  seg_mna_binned_tb %>%
  filter(seqnames == as.character(seqnames(view_granges))) %>% 
  group_by(start) %>%
  summarise(PercentSamples = 100*mean(isAmp),
            nSamples = sum(isAmp)) %>%
  ungroup() %>%
  ggplot(aes(x=start, y=PercentSamples)) +
  geom_line() + 
  ylab("Samples [%]") + 
  xlim(start(view_granges),end(view_granges)) +
  theme_kons2() 
print(real_profiles.fig)

ggsave("/Volumes/Elements/MYCNAmplicon/Results/MB_AggregateCopyNumberProfile.pdf",
       egg::ggarrange(
         genes.fig+
           theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
         #cRE.fig+
         #   theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
         real_profiles.fig,
         nrow = 2, 
         heights = c(0.05,1)),
       height=2, width=3, onefile = FALSE, useDingbats = F)


# ------------------------------------------------------------------------------
# Read ChIP-seq data
# ------------------------------------------------------------------------------

standardchrs = standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)[1:23]
chip_binsize = 1000 
chip_bins = tileGenome(seqlengths = seqlengths(BSgenome.Hsapiens.UCSC.hg19),tilewidth = chip_binsize, cut.last.tile.in.chrom = T) %>%
  filter_by_overlaps(view_granges, minoverlap = chip_binsize/10)
seqlevels(chip_bins, pruning.mode="coarse") = standardchrs

mb_chip_metadata = 
  data.table::fread("/Volumes/Elements/MYCNAmplicon/Data/LinEtAl_MB_ChIPseq/mb_chipseq_metadata.txt") %>% as_tibble() %>%
  filter(ChIP == "H3K27Ac")

# Read histone files
h3k27ac = list()
for (i in 1:nrow(mb_chip_metadata)){
  h3k27ac[[i]] = read_bigwig(mb_chip_metadata[[i,"File"]], overlap_ranges = view_granges)
  seqlevels(h3k27ac[[i]] , pruning.mode="coarse") = standardchrs
  h3k27ac[[i]] = binnedAverage(chip_bins, coverage(h3k27ac[[i]], weight=h3k27ac[[i]]$score), "h3k27ac")
  h3k27ac[[i]]$Sample = mb_chip_metadata[[i,"Sample"]]
  h3k27ac[[i]]$Subgroup = mb_chip_metadata[[i,"Subgroup"]]
}
h3k27ac = do.call(c,h3k27ac)

as_tibble(h3k27ac) %>% 
  mutate(Sample = as.character(Sample)) %>% 
  group_by(Subgroup) %>% 
  summarise(n=n_distinct(Sample))
# # A tibble: 4 x 2
# Subgroup     n
# <chr>    <int>
# 1 GROUP3      12
# 2 GROUP4      11
# 3 SHH          5
# 4 WNT          3

mean_h3k27ac = as_tibble(h3k27ac) %>%
  group_by(seqnames, start, end, Subgroup) %>%
  summarise(h3k27ac = mean(h3k27ac, na.rm=T)) %>%
  ungroup()

mean_h3k27ac %>%
  filter(Subgroup == "GROUP4") %>%
  mutate(value = as.numeric(h3k27ac)) %>% 
  dplyr::select(-Subgroup, -h3k27ac) %>% 
  write.table("/Volumes/Elements/MYCNAmplicon/Results/SourceData/MB_GROUP4_MeanH3K27ac_MYCN_neighborhood.txt",
              col.names=T, row.names=F, quote = F, sep = "\t")

mean_h3k27ac %>%
  filter(Subgroup == "SHH") %>%
  mutate(value = as.numeric(h3k27ac)) %>% 
  dplyr::select(-Subgroup, -h3k27ac) %>% 
  write.table("/Volumes/Elements/MYCNAmplicon/Results/SourceData/MB_SHH_MeanH3K27ac_MYCN_neighborhood.txt",
              col.names=T, row.names=F, quote = F, sep = "\t")

h3k27ac.fig = as_tibble(h3k27ac) %>%
  filter(Subgroup %in% c("GROUP4", "SHH")) %>% 
  group_by(start, Subgroup) %>%
  summarise(h3k27ac = mean(h3k27ac, na.rm=T)) %>%
  ungroup() %>% 
  ggplot(aes(x=start, y=h3k27ac, color=Subgroup)) + 
  geom_line() +
  xlim(start(view_granges), end(view_granges)) +
  facet_grid(Subgroup ~ .) + 
  guides(color=F) +
  scale_color_manual(values = c("SHH" = "#4daf4a", "GROUP4" = "#984ea3", "GROUP3"="#ff7f00", "Unknown"="grey70")) +
  theme_kons2() 
print(h3k27ac.fig)  

real_profiles_bygroup.fig =
  seg_mna_binned_tb %>%
  filter(Subgroup %in% c("GROUP4", "SHH")) %>% 
  filter(seqnames == as.character(seqnames(view_granges))) %>% 
  group_by(start, Subgroup) %>%
  summarise(PercentSamples = 100*mean(isAmp),
            nSamples = sum(isAmp)) %>%
  ungroup() %>%
  ggplot(aes(x=start, y=PercentSamples, color=Subgroup)) +
  geom_line() + 
  ylab("Samples [%]") + 
  xlim(start(view_granges),end(view_granges)) +
  scale_color_manual(values = c("SHH" = "#4daf4a", "GROUP4" = "#984ea3", "GROUP3"="#ff7f00", "Unknown"="grey70")) +
  guides(color=F) +
  theme_kons2() 
print(real_profiles_bygroup.fig)

ggsave("/Volumes/Elements/MYCNAmplicon/Results/MB_AggregateCopyNumberProfile_bySubgroup.pdf",
       egg::ggarrange(
         genes.fig
          +theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
         h3k27ac.fig 
          +theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
         real_profiles_bygroup.fig,
         nrow = 3, 
         heights = c(0.05,0.5,1)),
       height=3, width=3, onefile = FALSE, useDingbats = F)

seg_mna_tb %>%
  mutate(SubgroupSample = paste0(Subgroup, "_", Sample)) %>% 
  filter(seqnames == as.character(seqnames(view_granges))) %>% 
  ggplot(aes(x=start, y=SubgroupSample, color = Subgroup)) + 
  geom_errorbarh(aes(xmin=start, xmax=end), height=0) +
  #facet_grid(Subgroup ~ . ) +
  scale_color_manual(values = c("SHH" = "#4daf4a", "GROUP4" = "#984ea3", "GROUP3"="#ff7f00", "Unknown"="grey70")) +
  theme_kons2() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y=element_blank()) +
  ylab("") +
  xlab(seqnames(view_granges)) +
  annotate(geom = "rect", 
           xmin=start(mycn_granges),
           xmax=end(mycn_granges),
           ymin=0, ymax=Inf) +
  facet_zoom(xlim = c(start(view_granges), end(view_granges))) +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ggtitle("Medulloblastoma MYCN Amplicons") +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/mb-arrays_MYCN_amplicons_bySubgroup.pdf",
         height=2, width=3, onefile=F, useDingbats=F)

# with sample names ....
seg_mna_tb %>%
  mutate(SubgroupSample = paste0(Subgroup, "_", Sample)) %>% 
  filter(seqnames == as.character(seqnames(view_granges))) %>% 
  ggplot(aes(x=start, y=SubgroupSample, color = Subgroup)) + 
  geom_errorbarh(aes(xmin=start, xmax=end), height=0) +
  #facet_grid(Subgroup ~ . ) +
  scale_color_manual(values = c("SHH" = "#4daf4a", "GROUP4" = "#984ea3", "GROUP3"="#ff7f00", "Unknown"="grey70")) +
  theme_kons2() +
  #theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y=element_blank()) +
  ylab("") +
  xlab(seqnames(view_granges)) +
  annotate(geom = "rect", 
           xmin=start(mycn_granges),
           xmax=end(mycn_granges),
           ymin=0, ymax=Inf) +
  facet_zoom(xlim = c(start(view_granges), end(view_granges))) +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  ggtitle("Medulloblastoma MYCN Amplicons")

# The SHH-WNT super enhancer that explains the plateau in SHH is at chr2:16437375-16507828
# source: http://genome.ucsc.edu/cgi-bin/hgc?hgsid=812406981_aRyQERdxHP2nOHPM6pAIQAY5udaM&c=chr2&l=15600380&r=16567430&o=16437374&t=16507828&g=hub_74989_SHH%2DWNT_SE_bb&i=1_shh_859_lociStitched

# there are further SHH(-WNT)-SE surrounding MYCN:
# chr2:16152105-16211766

# ------------------------------------------------------------------------------
# one very short SHH: Sample MB-1179
# ------------------------------------------------------------------------------

# Look at 

# seg_mna_tb %>%
#   filter(Sample == "MB-1179") %>% View

# There is one very short amplicon that only overlaps a part of MYCN
# chr2	16052754	16083690

#> seg_mna_tb %>%
#  +        filter(Sample == "MB-1179") # that would also show non-chr2 amplicons
# A tibble: 8 x 8
# seqnames    start      end  width strand Value Sample  Subgroup
# <fct>       <int>    <int>  <int> <fct>  <dbl> <chr>   <chr>   
# 1 chr2     14337712 14526593 188882 *       2.40 MB-1179 SHH   # no Conserved SE, no SHH SE, but a enhancer in SHH and WNT samples 
# 2 chr2     14631081 14679093  48013 *       2.37 MB-1179 SHH   # no Conserved SE, no SHH SE, but clear enhancer in Group3 and Group4  
# 3 chr2     14777116 14788164  11049 *       2.77 MB-1179 SHH   # no Conserved SE, no SHH SE
# 4 chr2     14977816 15026488  48673 *       2.43 MB-1179 SHH   # no Conserved SE, no SHH SE   
# 5 chr2     15035998 15137809 101812 *       2.37 MB-1179 SHH   # no Conserved SE, no SHH SE  
# 6 chr2     15210238 15438727 228490 *       2.40 MB-1179 SHH   # no Conserved SE, no SHH SE
# 7 chr2     16052754 16083690  30937 *       2.31 MB-1179 SHH   # first coding exon of MYCN
# 8 chr2     16120135 16152517  32383 *       2.28 MB-1179 SHH   # no Conserved SE, no SHH SE, includes enhancer in many medulloblastoma samples across subgroups

# ------------------------------------------------------------------------------
# very short GROUP4 amplicons 
# ------------------------------------------------------------------------------

# one *very* short: MB-898
seg_mna_tb %>%
     filter(Sample == "MB-898") # would also show non-chr2 amplicons

# complex amplicon with 5 amplified segments on chr2
# A tibble: 5 x 8
# seqnames    start      end  width strand Value Sample Subgroup
# <fct>       <int>    <int>  <int> <fct>  <dbl> <chr>  <chr>   
# 1 chr2     15329510 15355850  26341 *       2.39 MB-898 GROUP4  # no Conserved SE, no GROUP4 SE, but an intronic H3K27ac peak in NBAS in some GROUP4 samples
# 2 chr2     16069542 16113652  44111 *       2.08 MB-898 GROUP4  # inludes MYCN, but no local enhancers
# 3 chr2     32786492 32819940  33449 *       2.40 MB-898 GROUP4  # no Conserved SE, no GROUP4 SE, includes parts of BIRC6 gene, small H3K27ac peak in few GROUP3 (sic) samples
# 4 chr2     58411336 58518273 106938 *       2.46 MB-898 GROUP4  # no Conserved SE, no GROUP4 SE, active FANCL
# 5 chr2     66406496 66451782  45287 *       2.19 MB-898 GROUP4  # no Conserved SE, no GROUP4 SE

# one that lacks most of the upstream part: sample MB-1238
# > seg_mna_tb %>%
#   +     filter(Sample == "MB-1238")
# # A tibble: 1 x 8
# seqnames    start      end  width strand Value Sample  Subgroup
# <fct>       <int>    <int>  <int> <fct>  <dbl> <chr>   <chr>   
#   1 chr2     16032278 16716234 683957 *       2.42 MB-1238 GROUP4  # ---->  includes all the downstream enhancers

# one that lacks most of the downstream part: sample MB-754
# > seg_mna_tb %>%
#   +     filter(Sample == "MB-754")
# # A tibble: 3 x 8
# seqnames    start      end  width strand Value Sample Subgroup
# <fct>       <int>    <int>  <int> <fct>  <dbl> <chr>  <chr>   
# 1 chr2     14827599 15163440 335842 *       2.30 MB-754 GROUP4  # no Conserved SE, no GROUP4 SE, but clear H3K27ac peak in GROUP4 samples and most other subgroups
# 2 chr2     15484066 15827604 343539 *       2.43 MB-754 GROUP4  # no Conserved SE, no GROUP4 SE, includes NBAS promotor and parts of gene and another clear H3K27ac peak in GROUP4 samples
# 3 chr2     15827604 16135550 307947 *       2.12 MB-754 GROUP4  # includes the upstream local enhancers of MYCN and the first downstream local enhancer

# third shortest upstream part : sample MB-1021
# > seg_mna_tb %>%
#   +     filter(Sample == "MB-1021")
# # A tibble: 1 x 8
# seqnames    start      end  width strand Value Sample  Subgroup
# <fct>       <int>    <int>  <int> <fct>  <dbl> <chr>   <chr>   
#   1 chr2     15978393 16304391 325999 *       1.80 MB-1021 GROUP4 # still includes the two local upstream enhancers

# ---> all but one (shortest, MB-898) sample include the first downstream enhancer 
# ---> all but two samples include the first upstream 
# ---> all but three samples include the second upstream


# # ------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------
# #                                   REST
# # ------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------
# # ------------------------------------------------------------------------------
# 
# 
# # ------------------------------------------------------------------------------
# # Now look samples with small MYCN-surronding amplicon parts
# # ------------------------------------------------------------------------------
# 
# for (i in 1:length(seg_mna)){
#   
#   this_seg_mna = seg_mna[[i]]
#   MYCN_frag_Value = seg_mna[[i]] %>% filter_by_overlaps(mycn_granges) %>% .$Value
#   MYCN_frag_start = seg_mna[[i]] %>% filter_by_overlaps(mycn_granges) %>% start()
#   MYCN_frag_end = seg_mna[[i]] %>% filter_by_overlaps(mycn_granges) %>% end()
#   
#   MYCN_frag_length = seg_mna[[i]] %>% filter_by_overlaps(mycn_granges) %>% width()
#   this_seg_mna = this_seg_mna %>% filter(abs(Value-MYCN_frag_Value)<0.3)
#   frags_n = this_seg_mna %>% length()
#   frags_lengthsum = this_seg_mna %>% width() %>% sum()
#   
#   if (i == 1){
#     amp_complexity = 
#       data.frame(
#         "Sample" = this_seg_mna$Sample[1],
#         "MYCNFragmentLength" = MYCN_frag_length,
#         "MYCNFragmentStart" = MYCN_frag_start,
#         "MYCNFragmentEnd" = MYCN_frag_end,
#         "NFrags" = frags_n,
#         "FragsLengths" = frags_lengthsum
#       )
#   } else {
#     amp_complexity = 
#       rbind(amp_complexity, 
#             data.frame(
#               "Sample" = this_seg_mna$Sample[1],
#               "MYCNFragmentLength" = MYCN_frag_length,
#               "MYCNFragmentStart" = MYCN_frag_start,
#               "MYCNFragmentEnd" = MYCN_frag_end,
#               "NFrags" = frags_n,
#               "FragsLengths" = frags_lengthsum
#             ))
#   }
# }
# 
# # Lets look at that
# amp_complexity %>% View
# # ---> there are two clear outliers in MYCN-surrounding fragment length --> both of them have additional amplicons
# 
# # Plot amplicon complexity vs. MYCN-surrounding fragment length
# amp_complexity %>%
#   ggplot(aes(x=MYCNFragmentLength, y=FragsLengths)) + 
#   geom_point() +
#   geom_abline() +
#   scale_x_log10() + 
#   scale_y_log10()
# 
# # Are complex amplicons associated with a short MYCN-surrounding fragment?
# amp_complexity = amp_complexity %>%
#   mutate(ComplexAmplicon = NFrags > 1)
# t.test(MYCNFragmentLength ~ ComplexAmplicon, data = amp_complexity) # p=0.08 nicht sign.
# 
# # 
# amp_complexity %>%
#   ggplot(aes(x=MYCNFragmentLength)) + 
#   geom_density()
# 
# amp_complexity %>%
#   ggplot(aes(x=FragsLengths)) + 
#   geom_density()
# mean(amp_complexity$FragsLengths)
# median(amp_complexity$FragsLengths)
# 
# # ------------------------------------------------------------------------------
# # investigate skew region
# # ------------------------------------------------------------------------------
# `%notin%` <- Negate(`%in%`)
# outlier_samples = c("GSM918858_MT792_MDT-MB-898", "GSM919158_MT1067_MDT-MB-1179")
# 
# seg_mna_tb = as_tibble(unlist(seg_mna))
# 
# # manually define local peak 
# seg_mna_binned_tb %>%
#   filter(seqnames == "chr2") %>% 
#   group_by(start) %>%
#   summarise(PercentSamples = 100*mean(isAmp),
#             nSamples = sum(isAmp)) %>%
#   ungroup() %>%
#   filter(PercentSamples > 40) %>%
#   ggplot(aes(x=start, y=PercentSamples)) +
#   geom_line()
# 
# seg_mna_binned_tb %>%
#   filter(seqnames == "chr2") %>% 
#   group_by(start) %>%
#   summarise(PercentSamples = 100*mean(isAmp),
#             nSamples = sum(isAmp)) %>%
#   ungroup() %>%
#   filter(PercentSamples > 40) %>%
#   View
# 
# skew_region = GRanges(
#   seqnames = "chr2",
#   ranges = IRanges(start = 16510001, end = 16520001)
# )
# skew_region = GRanges(
#   seqnames = "chr2",
#   ranges = IRanges(start = 15930001, end = 15940000)
# )
# skew_region = GRanges(
#   seqnames = "chr2",
#   ranges = IRanges(start = 16270001, end = 16280001)
# )
# 
# coamplify_skew = 
#   seg_mna_tb %>%
#   filter(hasOverlap_withChr(seqnames, start, end, as.character(seqnames(skew_region)), start(skew_region), end(skew_region))) %>%
#   .$Sample %>%
#   unique()
# length(coamplify_skew)
# 
# non_coamplify_skew = 
#   seg_mna_tb %>%
#   mutate(overlapwithskew = hasOverlap_withChr(seqnames, start, end, as.character(seqnames(skew_region)), start(skew_region), end(skew_region))) %>%
#   group_by(Sample) %>%
#   summarise(overlapwithskew = any(overlapwithskew)) %>%
#   ungroup() %>%
#   filter(!overlapwithskew) %>%
#   .$Sample %>%
#   unique()
# length(non_coamplify_skew)
# 
# # compare profiles for class I vs class II
# margin_of_interest = 1000000
# view_chr = as.character(seqnames(mycn_granges))
# view_start = start(mycn_granges)-margin_of_interest
# view_end = end(mycn_granges)+margin_of_interest
# view_granges = GRanges(seqnames = view_chr, ranges = IRanges(start = view_start, end = view_end))
# 
# seg_mna_binned_tb %>%
#   #filter(Name %notin% outlier_samples) %>% 
#   mutate(coamplifiesSkew = Name %in% coamplify_skew) %>% 
#   filter(seqnames == "chr2") %>% 
#   group_by(start, coamplifiesSkew) %>%
#   summarise(PercentSamples = 100*mean(isAmp),
#             nSamples = sum(isAmp)) %>%
#   ungroup() %>%
#   ggplot(aes(x=start, y=PercentSamples, color=coamplifiesSkew)) +
#   geom_line() + 
#   ylab("Samples [%]") + 
#   xlim(view_start, view_end) +
#   xlab("chr2") +
#   theme_kons2() +
#   guides(color=F)
# 
# # compare complexity for class I vs class II
# seg_mna_tb %>%
#   group_by(Sample) %>%
#   summarise(n=n_distinct(start)) %>% 
#   ungroup() %>%
#   mutate(isComplex = n>1,
#          coamplifiesSkew = Sample %in% coamplify_skew) %>%
#   ggplot(aes(x=isComplex, fill=coamplifiesSkew, group=coamplifiesSkew)) + 
#   geom_bar(position="dodge")
# 
# # do a Fisher test.
# seg_mna_tb %>%
#   group_by(Sample) %>%
#   summarise(n=n_distinct(start)) %>% 
#   ungroup() %>%
#   mutate(isComplex = n>1,
#          coamplifiesSkew = Sample %in% coamplify_skew) %>%
#   group_by(coamplifiesSkew, isComplex) %>%
#   summarise(n = n_distinct(Sample)) %>%
#   .$n %>% matrix(nrow=2) %>%
#   fisher.test()
# ## --> no relationship between complexity and class I vs. class II
# 
# # ------------------------------------------------------------------------------
# # 
# # ------------------------------------------------------------------------------
# 
# single_amp_fragments = 
#   seg_mna_tb %>%
#   mutate(targetamp = hasOverlap_withChr(seqnames, start, end, "chr2", start(mycn_granges), end(mycn_granges))) %>%
#   group_by(Sample) %>%
#   summarise(nFragments = n_distinct(start),
#             targetamp = any(targetamp)) %>%
#   ungroup() %>%
#   filter(nFragments == 1, targetamp == TRUE) %>%
#   .$Sample %>%
#   unique()
# length(single_amp_fragments)
# 
# this_seg = 
#   seg_mna_tb %>%
#   filter(Sample %in% single_amp_fragments)
# 
# n_ampl = length(unique(this_seg$Sample))
# 
# seg_mna_binned_tb %>%
#   mutate(isSingleAmpFragment = Name %in% single_amp_fragments) %>% 
#   filter(seqnames == "chr2") %>% 
#   group_by(start,isSingleAmpFragment) %>%
#   summarise(PercentSamples = 100*mean(isAmp),
#             nSamples = sum(isAmp)) %>%
#   ungroup() %>%
#   ggplot(aes(x=start, y=PercentSamples, color=isSingleAmpFragment)) +
#   geom_line() + 
#   ylab("Samples [%]") + 
#   xlim(start(mycn_granges)-roi_margin, end(mycn_granges)+roi_margin) +
#   xlab("chr2") +
#   theme_kons2() +
#   guides(color=F)

