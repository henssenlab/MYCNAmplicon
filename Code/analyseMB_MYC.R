library(data.table)
library(dplyr)
library(plyranges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)
library(ggforce)
source("/Volumes/Elements/MYCNAmplicon/Code/CustomThemes.R")
source("/Volumes/Elements/MYCNAmplicon/Code/hasOverlap.R")

# MYC
gene_of_interest = GRanges(
  seqnames = c("chr8"),
  ranges = IRanges(
    start = c(128747680),
    end = c(128753674)
  )
)
gene_of_interest$hgnc_symbol = "MYC"

seg =
  data.table::fread("/Volumes/Elements/MYCNAmplicon/Data/mb-arrays.rawcopy.segments.txt") %>%
  as_tibble() %>%
  dplyr::select(-Allelic.Imbalance) %>%
  filter(Value >= 1.8) %>%
  dplyr::mutate(ID = Sample) %>%
  makeGRangesListFromDataFrame(split.field = "ID", keep.extra.columns=T)
seqlevels(seg, pruning.mode="coarse") = standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)

arrayDataIsAmp = vector()
for (sample_idx in 1:length(seg)){
  if (length(findOverlaps(gene_of_interest, seg[sample_idx])) > 0) arrayDataIsAmp = c(arrayDataIsAmp, sample_idx)
}
seg_amp = seg[arrayDataIsAmp]

# how many cases are gene_of_interest-amplified?
length(seg_amp)

# compute coverage
seg_amp_cov = lapply(seg_amp, coverage)

# Bin Everything For Plotting
bins = tileGenome(seqlengths = seqlengths(BSgenome.Hsapiens.UCSC.hg19),
                  tilewidth = 10000,
                  cut.last.tile.in.chrom = T)
seqlevels(bins, pruning.mode = "coarse") = standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)
seg_amp_binned = lapply(seg_amp_cov, function (this_profile_cov) binnedAverage(bins, numvar = this_profile_cov, varname = "isAmp"))
for (i in 1:length(seg_amp_binned)){
  seg_amp_binned[[i]]$Name = names(seg_amp)[i]
}
seg_amp_binned_tb = do.call(rbind, lapply(seg_amp_binned, as_tibble))
seg_amp_binned_tb$isAmp = seg_amp_binned_tb$isAmp > 0 # allow only 0 or 1 as bin values = all bins with amplicon overlap --> 1

margin_of_interest = 2000000
seg_amp_binned_tb %>%
  filter(seqnames == as.character(seqnames(gene_of_interest))) %>% 
  group_by(start) %>%
  summarise(PercentSamples = 100*mean(isAmp),
            nSamples = sum(isAmp)) %>%
  ungroup() %>%
  ggplot(aes(x=start, y=PercentSamples)) +
  geom_line() + 
  ylab("Samples [%]") + 
  xlim(start(gene_of_interest)-margin_of_interest,end(gene_of_interest) + margin_of_interest) +
  xlab(as.character(seqnames(gene_of_interest))) +
  theme_kons2() 

# look at single-fragment cases ------------------------------------------------

seg_amp_tb = as_tibble(unlist(seg_amp))
single_amp_fragments = 
  seg_amp_tb %>%
  mutate(targetamp = hasOverlap_withChr(seqnames, start, end, as.character(seqnames(gene_of_interest)), start(gene_of_interest), end(gene_of_interest))) %>%
  group_by(Sample) %>%
  summarise(nFragments = n_distinct(start),
            targetamp = any(targetamp)) %>%
  ungroup() %>%
  filter(nFragments == 1, targetamp == TRUE) %>%
  .$Sample %>%
  unique()

length(single_amp_fragments)

# Plotting
margin_of_interest = 5000000
view_chr = as.character(seqnames(gene_of_interest))
view_start = start(gene_of_interest)-margin_of_interest
view_end = end(gene_of_interest)+margin_of_interest

real_profiles.fig = 
  seg_amp_binned_tb %>%
  mutate(isSingleAmpFragment = Name %in% single_amp_fragments) %>% 
  filter(seqnames == view_chr) %>% 
  group_by(start,isSingleAmpFragment) %>%
  summarise(PercentSamples = 100*mean(isAmp),
            nSamples = sum(isAmp)) %>%
  ungroup() %>%
  ggplot(aes(x=start, y=PercentSamples, color=isSingleAmpFragment)) +
  geom_line() + 
  ylab("Samples [%]") + 
  xlim(view_start, view_end) +
  xlab(as.character(seqnames(gene_of_interest))) +
  theme_kons2() +
  guides(color=F)

ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl", host="grch37.ensembl.org")
genes_plotting_df =  getBM(attributes=c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "gene_biotype", "ensembl_gene_id"),
                           filters= "chromosome_name", 
                           values= gsub("chr", "", view_chr), 
                           mart=ensembl)
genes_plotting_df$my_name = ifelse(genes_plotting_df$hgnc_symbol != "", 
                                   genes_plotting_df$hgnc_symbol, 
                                   paste0(genes_plotting_df$ensembl_gene_id))
genes_plotting_df = genes_plotting_df %>%
  filter(gene_biotype == "protein_coding") %>% 
  filter(hasOverlap_withChr(gsub("chr", "", view_chr), view_start, view_end,
                            chromosome_name, start_position, end_position))
genes.fig = genes_plotting_df %>%
  ggplot(aes(x=start_position, y=1)) +
  geom_rect(xmin=genes_plotting_df$start_position, xmax=genes_plotting_df$end_position, ymin=-Inf, ymax=0.1, color=NA, fill="black", alpha=0.5) + 
  ylim(0,2) +
  geom_text_repel(aes(x=start_position, y=0.1, label = hgnc_symbol), nudge_y=0.1, size=2, segment.size=0.1, min.segment.length = 0) +
  theme_kons2() + 
  theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank()) +
  xlim(view_start,view_end)

print(egg::ggarrange(
  genes.fig+
    theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
  real_profiles.fig,
  nrow = 2, 
  heights = c(0.5,1)))

ggsave("/Volumes/Elements/MYCNAmplicon/Results/TEST.pdf",
       egg::ggarrange(
         genes.fig+
           theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
         real_profiles.fig,
         nrow = 2, 
         heights = c(0.5,1)),
       height=3, width=3, onefile = FALSE, useDingbats = F)


# look at region of skew - posisitive vs. rest (is similar to class I vs class II) ---------------------------------

# manually define local peak 
seg_amp_binned_tb %>%
  filter(seqnames == view_chr) %>% 
  group_by(start) %>%
  summarise(PercentSamples = 100*mean(isAmp),
            nSamples = sum(isAmp)) %>%
  ungroup() %>%
  filter(PercentSamples > 60) %>%
  View

skew_region = GRanges(
  seqnames = "chr8",
  ranges = IRanges(start = 128320001, end = 128330000)
)

coamplify_skew = 
  seg_amp_tb %>%
  filter(hasOverlap_withChr(seqnames, start, end, as.character(seqnames(skew_region)), start(skew_region), end(skew_region))) %>%
  .$Sample %>%
  unique()
length(coamplify_skew)

non_coamplify_skew = 
  seg_amp_tb %>%
  mutate(overlapwithskew = hasOverlap_withChr(seqnames, start, end, as.character(seqnames(skew_region)), start(skew_region), end(skew_region))) %>%
  group_by(Sample) %>%
  summarise(overlapwithskew = any(overlapwithskew)) %>%
  ungroup() %>%
  filter(!overlapwithskew) %>%
  .$Sample %>%
  unique()
length(non_coamplify_skew)

# compare profiles for class I vs class II
seg_amp_binned_tb %>%
  mutate(coamplifiesSkew = Name %in% coamplify_skew) %>% 
  filter(seqnames == view_chr) %>% 
  group_by(start, coamplifiesSkew) %>%
  summarise(PercentSamples = 100*mean(isAmp),
            nSamples = sum(isAmp)) %>%
  ungroup() %>%
  ggplot(aes(x=start, y=PercentSamples, color=coamplifiesSkew)) +
  geom_line() + 
  ylab("Samples [%]") + 
  xlim(view_start, view_end) +
  xlab(as.character(seqnames(gene_of_interest))) +
  theme_kons2() +
  guides(color=F)

# compare complexity for class I vs class II
seg_amp_tb %>%
  group_by(Sample) %>%
  summarise(n=n_distinct(start)) %>% 
  ungroup() %>%
  mutate(isComplex = n>1,
         coamplifiesSkew = Sample %in% coamplify_skew) %>%
  ggplot(aes(x=isComplex, fill=coamplifiesSkew, group=coamplifiesSkew)) + 
  geom_bar(position="dodge")

# do a Fisher test.
seg_amp_tb %>%
  group_by(Sample) %>%
  summarise(n=n_distinct(start)) %>% 
  ungroup() %>%
  mutate(isComplex = n>1,
         coamplifiesSkew = Sample %in% coamplify_skew) %>%
  group_by(coamplifiesSkew, isComplex) %>%
  summarise(n = n_distinct(Sample))
cont.matrix = 
  matrix(c(1,4,6,19), nrow=2)
fisher.test(cont.matrix)
## ---> no relationship between complexity and class I vs. class II
