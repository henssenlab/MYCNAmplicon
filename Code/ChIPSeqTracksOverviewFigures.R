library(ggplot2)
library(dplyr)
library(data.table)
library(plyranges)
library(BSgenome.Hsapiens.UCSC.hg19)

source("/Volumes/Elements/MYCNAmplicon/Code/CustomThemes.R")

rnaseq_fname = "/Volumes/Elements/MYCNAmplicon/Workspace/boeva_rnaseq.sizeFactorNormalized.txt"
metadata_rnaseq_fname = "/Volumes/Elements/MYCNAmplicon/Data/Boeva_RNAseq_Metadata.csv"
metadata_chipseq_fname = "/Volumes/Elements/nb-cl-integr/metadata_h3k27.csv"

gene_of_interest = "MYCN"
gene_of_interest_ensembl = "ENSG00000134323"
gene_of_interest_chr = "chr2"
gene_of_interest_start = 16080683
gene_of_interest_end = 16087129
window_of_interest = 1500000
roi = GRanges(
  seqnames = gene_of_interest_chr,
  ranges = IRanges(start = gene_of_interest_start - window_of_interest,
                   end = gene_of_interest_end + window_of_interest),
  strand = "*"
)

standardchrs = standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)[1:23]
binsize = 1000 
bins = tileGenome(seqlengths = seqlengths(BSgenome.Hsapiens.UCSC.hg19),tilewidth = binsize, cut.last.tile.in.chrom = T) %>%
  filter_by_overlaps(roi, minoverlap = binsize/10)
seqlevels(bins, pruning.mode="coarse") = standardchrs

genes_df = 
  data.frame(
    "start" = c(12856998, 14772810, 15307032, 15731302, 16060521, 16080683, 16190549, 16730727, 17691851, 17720393, 17845079, 17935125, 17997763, 18059114, 18735989),
    "end" = c(12882860, 14790933, 15701454, 15771235, 16076139, 16087129, 16225923, 16847599, 17699706, 17838285, 17981509, 17966632, 17998368, 18542882, 18741959),
    "name" = c("TRIB2", "FAM84A", "NBAS", "DDX1", "MYCNUT", "MYCN", "GACAT3", "FAM49A", "RAD51AP2", "VSNL1", "SMC6", "GEN1", "MSGN1", "KCNS3", "RDH14"),
    "class" = rep("gene", 15)
  )

# ------------------------------------------------------------------------------
# load data
# ------------------------------------------------------------------------------

# samples
metadata_rnaseq = read.table(metadata_rnaseq_fname, header=T, sep=";", comment.char="#")
metadata_chipseq = read.table(metadata_chipseq_fname, header=T, sep=";", comment.char="#")
metadata_chipseq = metadata_chipseq %>% mutate(CellType = as.character(CellType),
                                               H3K27ac_bw = as.character(H3K27ac_bw),
                                               Input_bw = as.character(Input_bw),
                                               H3K27ac_peaks = as.character(H3K27ac_peaks)) %>%
  filter(CellType != "")

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

# Read histone files
controls = list()
h3k27ac = list()
h3k27ac_FC = list()
h3k27ac_peaks = list()
h3k27ac_peaks_unbinned = list()
for (i in 1:nrow(metadata_chipseq)){
  print(i)
  
  controls[[i]] = read_bigwig(metadata_chipseq[[i,"Input_bw"]], overlap_ranges = roi)
  seqlevels(controls[[i]] , pruning.mode="coarse") = standardchrs
  controls[[i]] = binnedAverage(bins, coverage(controls[[i]], weight=controls[[i]]$score), "input")
  controls[[i]]$CellType = metadata_chipseq[[i,"CellType"]]
  
  h3k27ac[[i]] = read_bigwig(metadata_chipseq[[i,"H3K27ac_bw"]], overlap_ranges = roi)
  seqlevels(h3k27ac[[i]] , pruning.mode="coarse") = standardchrs
  h3k27ac[[i]] = binnedAverage(bins, coverage(h3k27ac[[i]], weight=h3k27ac[[i]]$score), "h3k27ac")
  h3k27ac[[i]]$CellType = metadata_chipseq[[i,"CellType"]]
  
  h3k27ac_FC[[i]] = 
    join_overlap_inner(controls[[i]], h3k27ac[[i]], minoverlap=binsize) %>% 
    mutate(FC = (h3k27ac+0.01)/(input+0.01))
  
  h3k27ac_peaks_unbinned[[i]] = read_narrowpeaks(metadata_chipseq[[i,"H3K27ac_peaks"]], overlap_ranges = roi)
  seqlevels(h3k27ac_peaks_unbinned[[i]] , pruning.mode="coarse") = standardchrs
  h3k27ac_peaks_unbinned[[i]]$CellType = metadata_chipseq[[i,"CellType"]]
  
  h3k27ac_peaks[[i]] = read_narrowpeaks(metadata_chipseq[[i,"H3K27ac_peaks"]], overlap_ranges = roi)
  seqlevels(h3k27ac_peaks[[i]] , pruning.mode="coarse") = standardchrs
  h3k27ac_peaks[[i]] = binnedAverage(bins, coverage(h3k27ac_peaks[[i]]), "hasPeaks") %>% mutate(hasPeaks = (hasPeaks>0))
  h3k27ac_peaks[[i]]$CellType = metadata_chipseq[[i,"CellType"]]
}
controls = do.call(c,controls)
h3k27ac = do.call(c,h3k27ac)
h3k27ac_FC = do.call(c,h3k27ac_FC)
h3k27ac_peaks = do.call(c,h3k27ac_peaks)
h3k27ac_peaks_unbinned = do.call(c,h3k27ac_peaks_unbinned)

expression_overview = rnaseq %>% 
  filter(Study == "Boeva", CellType != "NCC")  %>%
  group_by(CellType, MYCNStatus) %>% 
  summarise(Expr = mean(Expr, na.rm=T)) %>% 
  mutate(hasExpr = Expr>100,
         ExprClass = ifelse(Expr < 100, "no MYCN", ifelse(MYCNStatus != "MNA", "low MYCN", "MYCN-amplified"))) %>%
  ungroup() %>% 
  mutate(CellType = forcats::fct_reorder(CellType, Expr, mean)) %>% 
  mutate(ExprClass = factor(ExprClass, levels=c("no MYCN", "low MYCN", "MYCN-amplified")))

non_mna_celltypes = expression_overview %>% filter(ExprClass != "MYCN-amplified") %>% .$CellType %>% as.character() %>% unique() 
mycnexp_non_mna_celltypes = expression_overview %>% filter(ExprClass == "low MYCN") %>% .$CellType %>% as.character() %>% unique() 

expression_overview = rnaseq %>% 
  filter(Study == "Boeva", CellType != "NCC")  %>%
  group_by(CellType, MYCNStatus) %>% 
  summarise(Expr = mean(Expr, na.rm=T)) %>% 
  mutate(hasExpr = Expr>100,
         ExprClass = ifelse(Expr < 100, "no MYCN", ifelse(MYCNStatus != "MNA", "low MYCN", "MYCN-amplified"))) %>%
  ungroup() %>% 
  mutate(CellType = forcats::fct_reorder(CellType, Expr, mean)) %>% 
  mutate(ExprClass = factor(ExprClass, levels=c("no MYCN", "low MYCN", "MYCN-amplified")))


# H3K27ac for all cell lines
h3k27ac.fig = 
  as_tibble(h3k27ac) %>%
  inner_join(expression_overview, by="CellType") %>%  
  mutate(CellType = forcats::fct_reorder(CellType, Expr, mean)) %>% 
  mutate(bin = paste0(seqnames, "_", start, "_", end)) %>%
  group_by(CellType, ExprClass, bin, seqnames, start, end) %>%
  summarise(h3k27ac = mean(h3k27ac)) %>%
  ungroup() %>%
  ggplot(aes(x=start, y=h3k27ac, color=ExprClass)) +
  geom_col() + 
  facet_grid(CellType ~ ., scales = "free_y") +
  theme_kons2() + 
  annotate(xmin=genes_df$start,xmax=genes_df$end, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) +
  annotate(xmin=cRE$start,xmax=cRE$end, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="red", alpha=0.2) +
  scale_color_manual(values=c("no MYCN"="steelblue", "MYCN-amplified"="black", "low MYCN"="firebrick3")) +
  theme(strip.text.y = element_text(angle = 0, size=8)) + 
  xlim(start(roi),end(roi)) + 
  theme(strip.text = element_text(size = 8, face = "plain")) +
  theme(panel.spacing.y = unit(6, units="pt")) +
  xlab("") + ylab("H3K7ac") +
  guides(color=F) 

h3k27ac.fig +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/BoevaData_H3K27ac_allNBCellLines.pdf", 
         height = 6, width = 3, useDingbats=F)

# genes.fig = genes_df %>%
#   ggplot(aes(x=start, y=1)) +
#   geom_rect(xmin=genes_df$start, xmax=genes_df$end, ymin=-Inf, ymax=Inf, color=NA, fill="black") + 
#   theme_kons2() + 
#   theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank()) +
#   xlim(start(roi),end(roi))
# cRE = read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_nMNA_MYCN_Enhancers.bed")
# cRE = as_tibble(cRE)
# cRE.fig = cRE %>%
#   ggplot(aes(x=start, y=1)) +
#   geom_rect(xmin=cRE$start, xmax=cRE$end, ymin=-Inf, ymax=Inf, color="firebrick3", fill="firebrick3") + 
#   theme_kons2() + 
#   theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank()) +
#   xlim(start(roi),end(roi))
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
#   xlim(start(roi),end(roi))
# 
# 
# ggsave("/Volumes/Elements/MYCNAmplicon/Results/BoevaData_H3K27ac_allNBCellLines.pdf", 
#        egg::ggarrange(
#          genes.fig +
#            theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
#          ggplot() + theme_void(),
#          cRE.fig +
#            theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
#          ggplot() + theme_void(),
#          h3k27ac.fig,
#          nrow=5,
#          heights = c(0.5,0.2,0.5,0.2,10)
#        ),
#        height = 10, width = 3, useDingbats=F)

