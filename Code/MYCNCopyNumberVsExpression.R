library(ggplot2)
library(dplyr)
library(data.table)
library(biomaRt)
library(ggrepel)
library(plyranges)
library(ggpubr)
library(BSgenome.Hsapiens.UCSC.hg19)
source("/Volumes/Elements/MYCNAmplicon/Code/CustomThemes.R")

# ------------------------------------------------------------------------------
# which gene
# ------------------------------------------------------------------------------

gene_of_interest = "MYCN"
gene_of_interest_ensembl = "ENSG00000134323"
gene_of_interest_chr = "chr2"
gene_of_interest_start = 16080683
gene_of_interest_end = 16087129
gene_of_interest_gr = GRanges(
  seqnames = gene_of_interest_chr,
  ranges = IRanges(start = gene_of_interest_start,
                   end = gene_of_interest_end),
  strand = "*"
)

# ------------------------------------------------------------------------------
# Read enhancer positions
# ------------------------------------------------------------------------------

mycn_enh = read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_nMNA_MYCN_Enhancers.bed")
e4 = mycn_enh %>% filter(name == "e4")

# ------------------------------------------------------------------------------
# Read Copy Number Data
# ------------------------------------------------------------------------------

metadata_chipseq = data.frame(
  CellType = c("GICAN", "SH-EP", "SK-N-AS", "GIMEN", "SK-N-SH", "NB69", "SJNB12",
               "SH-SY5Y", "SJNB1", "SK-N-FI", "CLB-GA", "NB-EBc1",
               "LAN1", "CLB-PE", "SK-N-DZ", "CLB-CAR", "CLB-MA",
               "IMR32", "CHP212", "SJNB8", "TR14", "SK-N-BE2-C",
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

metadata_chipseq$output_bed_fname = paste0("/Volumes/Elements/nb-cl-chipseq-qdnaseq/Boeva_", metadata_chipseq$CellType, "_Input.trimmed.bwa_hg19.rmdup.bam.qdnaseq.bed")

# cn = GRangesList()
# for (i in 1:nrow(metadata_chipseq)){
#   bed = fread(metadata_chipseq[i, "output_bed_fname"], skip = 1) %>% as_tibble()
#   colnames(bed) = c("chromosome", "start", "end", "interval", "copynumber", "strand")
#   bed_gr = makeGRangesFromDataFrame(bed, keep.extra.columns = T)
#   start(bed_gr) = start(bed_gr) + 1
#   strand(bed_gr) = "*"
#   seqlevelsStyle(bed_gr) <- "UCSC"
#   seqlevels(bed_gr, pruning.mode = "coarse") = seqlevels(BSgenome.Hsapiens.UCSC.hg19)
#   seqlengths(bed_gr) = seqlengths(BSgenome.Hsapiens.UCSC.hg19)
#   bed_gr_gaps = gaps(bed_gr)
#   bed_gr_gaps = bed_gr_gaps %>% filter(strand == "*")
#   bed_gr_gaps$copynumber = NA
#   bed_gr = sort(c(bed_gr, bed_gr_gaps))
#   bed_gr = coverage(bed_gr, weight = "copynumber") %>% GRanges()
#   bed_gr$Sample = as.character(metadata_chipseq[i, "CellType"])
#   cn[[i]] = bed_gr
# }
# cn = do.call(c, cn)
# save.image("/Volumes/Elements/MYCNAmplicon/Results/boeva_input_cn.Rdata")

load("/Volumes/Elements/MYCNAmplicon/Results/boeva_input_cn.Rdata")
cn$Sample = gsub("-", "", cn$Sample)

mycn_cn = 
  cn %>%
  join_overlap_intersect(gene_of_interest_gr) %>% 
  as_tibble() %>%
  mutate(score = score * width/width(gene_of_interest_gr)) %>% 
  group_by(Sample) %>% 
  summarise(score = sum(score))

e4_coamplified = 
  cn %>%
  join_overlap_intersect(e4) %>% 
  as_tibble() %>% 
  filter(score.x > 5) %>% 
  .$Sample
  

# ------------------------------------------------------------------------------
# Read RNA-seq Data
# ------------------------------------------------------------------------------

metadata_rnaseq_fname = "/Volumes/Elements/MYCNAmplicon/Data/Boeva_RNAseq_Metadata.csv"
rnaseq_fname = "/Volumes/Elements/MYCNAmplicon/Workspace/boeva_rnaseq.sizeFactorNormalized.txt"
rnaseq = fread(rnaseq_fname) %>% as_tibble() # gives a warning, this is fine 
colnames(rnaseq)[1] = "ENSEMBLGene"
colnames(rnaseq) = gsub("Boeva_","", colnames(rnaseq))

rnaseq = rnaseq %>% 
  filter(ENSEMBLGene == gene_of_interest_ensembl) %>%
  dplyr::select(-ENSEMBLGene) %>% as.matrix() %>% t()
colnames(rnaseq) = "Expr"
rnaseq = data.frame(rnaseq) %>% 
  mutate(Sample = rownames(rnaseq)) %>%
  as_tibble() 

# ------------------------------------------------------------------------------
# Analyse MYCN Copy Number vs. Expression
# ------------------------------------------------------------------------------

cn_vs_rna = 
  inner_join(mycn_cn, rnaseq)

cn_vs_rna = 
  cn_vs_rna %>% 
  filter(score > 5) %>% # only MNA
  mutate(Class = ifelse(Sample %in% e4_coamplified, "Class I", "Class II")) %>% 
  mutate(ExprByCN = Expr / score)

cn_vs_rna %>% 
  filter(score > 5) %>% # only MNA
  ggplot(aes(x=score, y=Expr)) + 
  geom_smooth(method = "lm", color = "steelblue", fill = NA, size = 0.5) +
  geom_point() +
  stat_cor(method = "pearson", label.x.npc = "right", label.y.npc = "bottom", hjust=0.925, vjust = 2.5, size=4) +
  theme_kons1() + 
  xlab("MYCN Genomic Copy Number") + 
  ylab("MYCN Expression") + 
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/MYCN_CN_vs_Expr.pdf",
         height = 3, width = 3, useDingbats = F)

cn_vs_rna %>% 
  filter(score > 5) %>% # only MNA
  mutate(Class = ifelse(Sample %in% e4_coamplified, "Class I", "Class II")) %>% 
  mutate(ExprByCN = Expr / score) %>% 
  ggplot(aes(x=Class, y=ExprByCN)) + 
  geom_jitter(width = 0.1) +
  geom_text_repel(aes(label = ifelse(Class == "Class II", Sample, NA))) +
  theme_kons1() +
  xlab("") +
  ylab("Copy Number-Normalized\nMYCN Expression") +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/MYCN_CNNormalizedExpression_vs_Class.pdf",
         height = 3, width = 3, useDingbats = F)

cn_vs_rna %>% 
  filter(score > 5) %>% # only MNA
  mutate(Class = ifelse(Sample %in% e4_coamplified, "Class I", "Class II")) %>% 
  ggplot(aes(x=Class, y=score)) + 
  geom_jitter(width = 0.1) +
  geom_text_repel(aes(label = ifelse(Class == "Class II", Sample, NA))) +
  theme_kons1() +
  xlab("") +
  ylab("MYCN Copy Number Ratio") +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/MYCN_CN_vs_Class.pdf",
         height = 3, width = 3, useDingbats = F)

cn_vs_rna %>% 
  filter(score > 5) %>% # only MNA
  mutate(Class = ifelse(Sample %in% e4_coamplified, "Class I", "Class II")) %>% 
  ggplot(aes(x=score, y=Expr, color = Class)) + 
  geom_smooth(method="lm") +
  geom_point() +
  stat_cor(method = "pearson", label.x.npc = "right", label.y.npc = "bottom", hjust=0.925, vjust = 2.5, size=4) +
  theme_kons1() + 
  xlab("MYCN Genomic Copy Number") + 
  ylab("MYCN Expression") + 
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/MYCN_CN_vs_Expr_vs_Class.pdf",
         height = 3, width = 3, useDingbats = F)

#wilcox_test(Expr ~ as.factor(Class), data = cn_vs_rna, distribution = "exact")
#wilcox_test(score ~ as.factor(Class), data = cn_vs_rna, distribution = "exact")


wilcox_test(ExprByCN ~ as.factor(Class), data = cn_vs_rna, distribution = "exact")

library(rankFD)
rank.two.samples(ExprByCN ~ as.factor(Class), data = cn_vs_rna)
