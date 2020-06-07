library(ggplot2)
library(dplyr)
library(data.table)
library(biomaRt)
library(ggrepel)

source("/Volumes/Elements/MYCNAmplicon/Code/CustomThemes.R")

rnaseq_fname = "/Volumes/Elements/MYCNAmplicon/Workspace/boeva_rnaseq.sizeFactorNormalized.txt"
metadata_rnaseq_fname = "/Volumes/Elements/MYCNAmplicon/Data/Boeva_RNAseq_Metadata.csv"
metadata_chipseq_fname = "/Volumes/Elements/nb-cl-integr/metadata_h3k27.csv"

gene_of_interest = "ALK"
gene_of_interest_ensembl = "ENSG00000171094"

# ------------------------------------------------------------------------------
# load data
# ------------------------------------------------------------------------------

# samples
metadata_rnaseq = read.table(metadata_rnaseq_fname, header=T, sep=";", comment.char="#")

# RNA-seq
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
rnaseq = fread(rnaseq_fname) %>% as_tibble() # gives a warning, this is fine 
colnames(rnaseq)[1] = "ENSEMBLGene"
rnaseq = rnaseq %>% 
  filter(ENSEMBLGene == gene_of_interest_ensembl) %>%
  dplyr::select(-ENSEMBLGene) %>% as.matrix() %>% t()
colnames(rnaseq) = "Expr"
rnaseq = data.frame(rnaseq) %>% 
  mutate(Sample = rownames(rnaseq)) %>%
  as_tibble() %>% 
  full_join(metadata_rnaseq) %>% 
  dplyr::select(-featureCounts_fname)

rnaseq %>% 
  filter(Study == "Boeva", CellType != "NCC")  %>%
  group_by(CellType, MYCNStatus) %>% 
  summarise(Expr = mean(Expr, na.rm=T)) %>% 
  ungroup() %>% 
  mutate(CellType = forcats::fct_reorder(CellType, Expr, mean)) %>% 
  ggplot(aes(x=CellType, y=Expr, color = (MYCNStatus == "MNA"))) + 
  geom_point() +
  scale_color_manual(values=c("TRUE" = "firebrick3", "FALSE"="steelblue")) + 
  guides(color=F) +
  theme_kons1() +
  scale_y_log10() + 
  xlab("") +
  ylab("ALK Expression") +
  theme(legend.title=element_blank()) +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/BoevaData_ALK_Expression.pdf", 
         height = 3, width = 8,
         useDingbats=F)

set.seed(58)
rnaseq %>% 
  filter(Study == "Boeva", CellType != "NCC")  %>%
  group_by(CellType, MYCNStatus) %>% 
  summarise(Expr = mean(Expr, na.rm=T)) %>% 
  ungroup() %>% 
  mutate(MYCNStatus = ifelse(MYCNStatus == "MNA", "MYCN-amplified", "not MYCN-amplified")) %>% 
  mutate(CellType = forcats::fct_reorder(CellType, Expr, mean)) %>% 
  ggplot(aes(x=MYCNStatus, y=Expr, color = (CellType == "IMR32"))) + 
  geom_jitter(width=0.15) +
  scale_color_manual(values=c("TRUE" = "firebrick3", "FALSE"="grey")) + 
  geom_text_repel(aes(label = ifelse(CellType == "IMR32", "IMR32", NA)), nudge_x=0.4, nudge_y=.1) +
  guides(color=F) +
  theme_kons1() +
  scale_y_log10() + 
  xlab("") +
  ylab("ALK Expression \n[sizeFactor-normalized read counts]") +
  theme(legend.title=element_blank()) +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/BoevaData_ALK_Expression_byMNA.pdf", 
         height = 3, width = 3,
         useDingbats=F)
 
