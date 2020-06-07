rm(list=ls())

library(dplyr)
library(data.table)
library(ggplot2)
library(plyranges)
library(biomaRt)
library(ggrepel)
source("/Volumes/Elements/MYCNAmplicon/Code/CustomThemes.R")
source("/Volumes/Elements/MYCNAmplicon/Code/hasOverlap.R")


angiosarc = fread("/Volumes/Elements/MYCNAmplicon/Data/MannerEtAl_Angiosarcoma/Manner_Angiosarcoma_SupplementalTable_Reformatted.csv", dec = ",") %>% as_tibble()

# how many samples overall?
angiosarc %>%
  dplyr::select(-CloneMidpoint, -CloneName, -Chromosome, -Mapping, -ChromosomeCen, -GeneIDs, -GeneSymbols, -GeneNames) %>% 
  ncol()

angiosarc %>%
  filter(Chromosome == "8") %>%
  dplyr::select(-CloneName, -Chromosome, -Mapping, -ChromosomeCen, -GeneIDs, -GeneSymbols, -GeneNames) %>%
  tidyr::gather("Sample", "Value", 2:22) %>%
  ggplot(aes(x = CloneMidpoint, y=Value, color=Sample)) + 
  geom_line()
# ---> no far-away amplicons on chr8, can plot locally only

# choose threshold to match number of 8 MYC-amplified samples as reported in paper
myc_amplified_samples = angiosarc %>%
  filter(GeneSymbols == "MYC;AC103819.3;AC084123.9;AC103819.3") %>%
  dplyr::select(-CloneMidpoint, -CloneName, -Chromosome, -Mapping, -ChromosomeCen, -GeneIDs, -GeneSymbols, -GeneNames) %>% 
  tidyr::gather("Sample", "Value", 1:21) %>%
  filter(Value > 1.4) %>%
  .$Sample %>%
  unique()
myc_amplified_samples

angiosarc %>%
  filter(Chromosome == "8") %>%
  dplyr::select(-CloneName, -Chromosome, -Mapping, -ChromosomeCen, -GeneIDs, -GeneSymbols, -GeneNames) %>%
  tidyr::gather("Sample", "Value", 2:22) %>%
  filter(Sample %in% myc_amplified_samples) %>% 
  filter(CloneMidpoint>125000000, CloneMidpoint < 132500000) %>% 
  ggplot(aes(x = CloneMidpoint, y=Value, color=Sample)) + 
  geom_point(size=0.5) + 
  geom_vline(xintercept=128810204, linetype = "dashed") + #MYC-containing clone
  guides(color=F) + 
  facet_grid(Sample ~ .) + 
  scale_color_brewer(palette = "Set2") +
  theme_kons2() +
  theme(strip.text.y = element_text(angle = 0)) +
  scale_y_continuous(breaks = c(0,4)) +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/MannerEtAl_MYCAmplifiedSamples.pdf",
         height=3, width=3, useDingbats=F)


angiosarc %>%
  filter(Chromosome == "8") %>%
  dplyr::select(-CloneName, -Chromosome, -Mapping, -ChromosomeCen, -GeneIDs, -GeneSymbols, -GeneNames) %>%
  tidyr::gather("Sample", "Value", 2:22) %>%
  filter(Sample %in% myc_amplified_samples) %>% 
  filter(CloneMidpoint>125000000, CloneMidpoint < 132500000) %>% 
  mutate(CloneMidpoint = as.factor(CloneMidpoint)) %>% 
  filter(Value > 1.4) %>%
  group_by(CloneMidpoint) %>%
  summarise(n = n_distinct(Sample)) %>%
  tidyr::complete(CloneMidpoint, fill=list(n=0)) %>%
  mutate(PercentAmplified=100*n/8) %>%
  mutate(CloneMidpoint = as.numeric(as.character(CloneMidpoint))) %>%
  View

# how large is the co-amplification region?
# 128747680 (Start MYC) - 128192436 (Midpoint furthest clone that is 100% co-amplified) = 555244

# Plot the genes within the region of interest
roi_gr = GRanges(
  seqnames = "8",
  ranges = IRanges(
    start = 125000000,
    end = 132500000
  )
)
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl", host="grch37.ensembl.org")
genes_plotting_df =  getBM(attributes=c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "gene_biotype", "ensembl_gene_id"),
                           filters= "chromosome_name", 
                           values= gsub("chr", "", seqnames(roi_gr)), 
                           mart=ensembl)
genes_plotting_df$my_name = ifelse(genes_plotting_df$hgnc_symbol != "", 
                                   genes_plotting_df$hgnc_symbol, 
                                   paste0(genes_plotting_df$ensembl_gene_id))
genes_plotting_df$hgnc_symbol = 
  ifelse(genes_plotting_df$hgnc_symbol %in% c( 
                                              "FAM84B",
                                              "POU5F1B", 
                                              "MYC", 
                                              "TMEM75"), 
         genes_plotting_df$hgnc_symbol, 
         NA)

genes_plotting_df = genes_plotting_df %>%
  filter(gene_biotype == "protein_coding") %>% 
  filter(hasOverlap_withChr(gsub("chr", "", as.character(seqnames(roi_gr))), start(roi_gr), end(roi_gr),
                            chromosome_name, start_position, end_position))
genes.fig = genes_plotting_df %>%
  ggplot(aes(x=start_position, y=1)) +
  geom_rect(xmin=genes_plotting_df$start_position, xmax=genes_plotting_df$end_position, ymin=-Inf, ymax=0.1, color=NA, fill="black", alpha=0.5) + 
  ylim(0,2) +
  geom_text_repel(aes(x=start_position, y=0.1, label = hgnc_symbol), nudge_y=0.1, size=2, segment.size=0.1, min.segment.length = 0) +
  theme_kons2() + 
  theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank()) +
  xlim(start(roi_gr),end(roi_gr))
print(genes.fig)

f = angiosarc %>%
  filter(Chromosome == "8") %>%
  dplyr::select(-CloneName, -Chromosome, -Mapping, -ChromosomeCen, -GeneIDs, -GeneSymbols, -GeneNames) %>%
  tidyr::gather("Sample", "Value", 2:22) %>%
  filter(Sample %in% myc_amplified_samples) %>% 
  filter(CloneMidpoint>125000000, CloneMidpoint < 132500000) %>% 
  mutate(CloneMidpoint = as.factor(CloneMidpoint)) %>% 
  filter(Value > 1.4) %>%
  group_by(CloneMidpoint) %>%
  summarise(n = n_distinct(Sample)) %>%
  tidyr::complete(CloneMidpoint, fill=list(n=0)) %>%
  mutate(PercentAmplified=100*n/8) %>%
  mutate(CloneMidpoint = as.numeric(as.character(CloneMidpoint))) %>% 
  ggplot(aes(x = CloneMidpoint, y=PercentAmplified)) + 
  geom_line(size=0.5) + 
  geom_vline(xintercept=128810204, linetype = "dashed") + #MYC-containing clone
  theme_kons2() +
  xlab("chromosome 8") + 
  ylab("Samples [%]") +
  xlim(start(roi_gr),end(roi_gr))

ggsave(
  "/Volumes/Elements/MYCNAmplicon/Results/MannerEtAl_AggregateAmplificationProfile.pdf",
  egg::ggarrange(
    genes.fig +
      ggtitle("Manner et al. Angiosarcoma (n=9)") +
      theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
    f,
    nrow=2,
    heights=c(0.5,2)
  ),
  height = 3, width = 3, useDingbats = F, onefile = F
)



#+
#  ggsave("/Volumes/Elements/MYCNAmplicon/Results/MannerEtAl_AggregateAmplificationProfile.pdf",
#         height=2, width=2, useDingbats=F)


