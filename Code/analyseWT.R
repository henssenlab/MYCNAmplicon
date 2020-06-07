rm(list=ls())

library(data.table)
library(dplyr)
library(plyranges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)
library(ggforce)
source("/Volumes/Elements/MYCNAmplicon/Code/CustomThemes.R")

# there are 125 donors but 130 samples in this data
seg = 
  data.table::fread("/Volumes/Elements/TARGET-WT/copy_number_somatic_mutation.WT-US.tsv") %>%
  as_tibble() %>%
  dplyr::select(icgc_donor_id, icgc_sample_id, segment_mean,	chromosome,	chromosome_start,	chromosome_end)
length(unique(seg$icgc_donor_id))

seg = 
  data.table::fread("/Volumes/Elements/TARGET-WT/copy_number_somatic_mutation.WT-US.tsv") %>%
  as_tibble() %>%
  dplyr::select(icgc_donor_id, icgc_sample_id, segment_mean,	chromosome,	chromosome_start,	chromosome_end) %>% 
  dplyr::rename(start=chromosome_start, end=chromosome_end, value = segment_mean) %>%
  mutate(chromosome = paste0("chr", chromosome)) %>% 
  filter(value >= 1.8) %>% 
  makeGRangesListFromDataFrame(split.field = "icgc_sample_id", keep.extra.columns=T)
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
  write_bed("/Volumes/Elements/MYCNAmplicon/Results/WT_AllCopyNumberProfiles.bed")
seg_mna %>%
  write_bed("/Volumes/Elements/MYCNAmplicon/Results/WT_MNANumberProfiles.bed")

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

# ------------------------------------------------------------------------------
# Figure
# ------------------------------------------------------------------------------
gene_of_interest = "MYCN"
gene_of_interest_ensembl = "ENSG00000134323"
gene_of_interest_chr = "chr2"
gene_of_interest_start = 16080683
gene_of_interest_end = 16087129
window_of_interest = 1250000

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
  xlim(gene_of_interest_start-1250000,gene_of_interest_end + 1250000)

# cRE = read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_nMNA_MYCN_Enhancers.bed") %>% as_tibble()
# cRE.fig = cRE %>%
#   ggplot(aes(x=start, y=score)) +
#   geom_rect(aes(xmin=start, xmax=end), ymin=-Inf, ymax=Inf, color="firebrick3", fill="firebrick3") + 
#   theme_kons2() + 
#   theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank()) +
#   xlim(gene_of_interest_start-500000,gene_of_interest_end + 500000)
# cRE.fig 
# cRE$class = "cRE"
# cRE_genes_plot_df = 
#   bind_rows(cRE, genes_df) %>% 
#   mutate(dummy=1) %>% 
#   dplyr::select(start, end, class, dummy) %>%
#   mutate(this_color = ifelse(class == "gene", NA, "firebrick3"),
#          this_fill = ifelse(class == "gene", "black", "firebrick3"))
# cRE_genes.fig = 
#   cRE_genes_plot_df %>%
#   ggplot(aes(x=start, y=dummy)) +
#   geom_rect(xmin=cRE_genes_plot_df$start, xmax=cRE_genes_plot_df$end, ymin=-Inf, ymax=Inf, color=cRE_genes_plot_df$this_color, fill=cRE_genes_plot_df$this_fill) + 
#   theme_kons2() +
#   theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank()) +
#   xlim(gene_of_interest_start-500000,gene_of_interest_end + 500000)

real_profiles.fig =
  seg_mna_binned_tb %>%
  filter(seqnames == "chr2") %>% 
  group_by(start) %>%
  summarise(PercentSamples = 100*mean(isAmp),
            nSamples = sum(isAmp)) %>%
  ungroup() %>%
  ggplot(aes(x=start, y=PercentSamples)) +
  geom_line() + 
  ylab("Samples [%]") + 
  xlim(gene_of_interest_start-1250000,gene_of_interest_end + 1250000) +
  theme_kons2() 

seg_mna_binned_tb %>%
  .$Name %>%
  unique() %>%
  length()
# ----> 16 MNA-WT samples  
  
# randomized_profiles.fig =
#   randomized_profiles_binned_tb %>%
#   group_by(start) %>%
#   summarise(PercentSamples = 100*mean(isAmp),
#             nSamples = sum(isAmp)) %>%
#   ungroup() %>%
#   ggplot(aes(x=start, y=PercentSamples)) +
#   geom_line(linetype = "dashed") + 
#   ylab("Samples [%]") + 
#   xlim(gene_of_interest_start-500000,gene_of_interest_end + 500000) +
#   theme_kons2() 

# mna_profiles_binned_tb$RealOrRandomized = "Real"
# randomized_profiles_binned_tb$RealOrRandomized = "Randomized"
# real_and_randomized_profiles.fig = bind_rows(mna_profiles_binned_tb, randomized_profiles_binned_tb) %>%
#   group_by(start, RealOrRandomized) %>%
#   summarise(PercentSamples = 100*mean(isAmp),
#             nSamples = sum(isAmp)) %>%
#   ungroup() %>%
#   ggplot(aes(x=start, y=PercentSamples, linetype = RealOrRandomized)) +
#   geom_line() + 
#   theme_kons1() + 
#   theme(strip.background = element_rect(fill="grey90", color=NA)) +
#   ylab("Samples [%]") + 
#   scale_linetype_manual(values = c("Real"="solid", "Randomized"="dotted")) +
#   guides(linetype = F) +
#   scale_y_continuous(breaks = c(40, 60, 80, 100), limits=c(40,100), expand=c(0,0)) +
#   xlim(gene_of_interest_start-500000,gene_of_interest_end + 500000) +
#   xlab("chr2") +
#   theme_kons2() 

ggsave("/Volumes/Elements/MYCNAmplicon/Results/WT_AggregateCopyNumberProfile.pdf",
       egg::ggarrange(
         genes.fig+
           theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
         #cRE.fig+
         #   theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
         real_profiles.fig,
         nrow = 2, 
         heights = c(0.05,1)),
       height=2, width=3, onefile = FALSE, useDingbats = F)

## Is WT class II frequency different from class II frequency in neuroblastoma?
cont.matrix = matrix(c(24,216,0,16), nrow=2)
fisher.test(cont.matrix)


