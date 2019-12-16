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
window_of_interest = 500000
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

# Generate a table for that
deltaFC = 
  as_tibble(h3k27ac_FC) %>%
  mutate(CellType = CellType.x) %>% dplyr::select(-CellType.y, -CellType.x) %>% 
  inner_join(expression_overview, by="CellType") %>%  
  mutate(CellType = forcats::fct_reorder(CellType, Expr, mean)) %>% 
  mutate(bin = paste0(seqnames, "_", start, "_", end)) %>%
  group_by(ExprClass, bin, seqnames, start, end) %>%
  summarise(FC = mean(FC)) %>%
  ungroup() %>%
  filter(ExprClass != "MYCN-amplified") %>%
  tidyr::spread(ExprClass, FC) %>% 
  mutate(DeltaFC = `low MYCN` - `no MYCN`) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)

deltaFC_bigwig = deltaFC[,"DeltaFC"] %>%
 mutate(score = DeltaFC) 
seqlevels(deltaFC_bigwig, pruning.mode="coarse") = standardchrs
seqlengths(deltaFC_bigwig) = seqlengths(BSgenome.Hsapiens.UCSC.hg19)[standardchrs]

deltaFC_bigwig %>%
 write_bigwig("/Volumes/Elements/MYCNAmplicon/Results/Boeva_nMNA_deltaFC.bigwig")

peak_regions = 
  coverage(h3k27ac_peaks_unbinned %>% 
             filter(CellType %in% mycnexp_non_mna_celltypes)) %>%
  GRanges() %>%
  filter(score >= 3) %>% 
  reduce(min.gapwidth = 2000)

peak_regions_dFC = binnedAverage(bins = peak_regions, 
                                 numvar = coverage(deltaFC, weight = "DeltaFC"),
                                 varname = "meanDeltaFC")
peak_regions_dFC$maxDFC = 
  group_by_overlaps(peak_regions, deltaFC) %>%
  summarise(maxDFC = max(DeltaFC)) %>% 
  .$maxDFC

peak_regions_dFC %>% 
  as_tibble() %>%
  mutate(maxDFC_rank = rank(-maxDFC)) %>%
  ggplot(aes(x=maxDFC_rank, y = maxDFC, color = maxDFC_rank<=6)) +
  geom_point() + 
  theme_kons1() + 
  scale_color_manual(values = c("black", "firebrick3")) + 
  guides(color=F) + 
  xlab("") + ylab("Maximum Fold Change\nDifference") + 
  theme(axis.ticks.x=element_blank(), axis.text.x = element_blank()) +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/RankCREByMaxDeltaFC.pdf",
         height=3, width=3)

cRE = peak_regions_dFC %>% 
  as_tibble() %>%
  mutate(maxDFC_rank = rank(-maxDFC)) %>% 
  filter(maxDFC_rank <= 6) # entspricht maxDFC >= 7
cRE$ID = c("e1", "MYCNp", "e2", "e3", "e4", "e5")

cRE %>%
  dplyr::select(seqnames, start, end, ID, maxDFC) %>% 
  write.table(file = "/Volumes/Elements/MYCNAmplicon/Results/Boeva_nMNA_MYCN_Enhancers.bed",
              row.names = F, col.names = F, sep = "\t", quote = F)

rnaseq %>% 
  filter(Study == "Boeva", CellType != "NCC")  %>%
  group_by(CellType, MYCNStatus) %>% 
  summarise(Expr = mean(Expr, na.rm=T)) %>% 
  mutate(hasExpr = Expr>100,
         ExprClass = ifelse(Expr < 100, "no MYCN", ifelse(MYCNStatus != "MNA", "low MYCN", "MYCN-amplified"))) %>%
  ungroup() %>% 
  mutate(CellType = forcats::fct_reorder(CellType, Expr, mean)) %>% 
  mutate(ExprClass = factor(ExprClass, levels=c("no MYCN", "low MYCN", "MYCN-amplified"))) %>% 
  ggplot(aes(x=CellType, y=Expr, color=ExprClass)) + 
  geom_point() +
  scale_color_manual(values=c("no MYCN"="steelblue", "MYCN-amplified"="gold2", "low MYCN"="firebrick3")) +
  theme_kons1() +
  scale_y_log10() + 
  xlab("") +
  ylab("MYCN Expression") +
  theme(legend.title=element_blank()) +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/BoevaData_Expression.pdf", 
         height = 3, width = 8,
         useDingbats=F)

rnaseq %>% 
  filter(Study == "Boeva", CellType != "NCC")  %>%
  group_by(CellType, MYCNStatus) %>% 
  summarise(Expr = mean(Expr, na.rm=T)) %>% 
  mutate(hasExpr = Expr>100,
         ExprClass = ifelse(Expr < 100, "no MYCN", ifelse(MYCNStatus != "MNA", "low MYCN", "MYCN-amplified"))) %>%
  ungroup() %>% 
  mutate(CellType = forcats::fct_reorder(CellType, Expr, mean)) %>% 
  mutate(ExprClass = factor(ExprClass, levels=c("no MYCN", "low MYCN", "MYCN-amplified"))) %>% 
  filter(ExprClass != "MYCN-amplified") %>% 
  ggplot(aes(x=CellType, y=Expr, color=ExprClass)) + 
  geom_point() +
  scale_color_manual(values=c("no MYCN"="steelblue", "MYCN-amplified"="gold2", "low MYCN"="firebrick3")) +
  theme_kons1() +
  xlab("") +
  ylab("MYCN Expression") +
  theme(legend.title=element_blank()) +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/BoevaData_Expression_onlynonMNA.pdf", 
         height = 3, width = 4,
         useDingbats=F)

expression_overview = rnaseq %>% 
  filter(Study == "Boeva", CellType != "NCC")  %>%
  group_by(CellType, MYCNStatus) %>% 
  summarise(Expr = mean(Expr, na.rm=T)) %>% 
  mutate(hasExpr = Expr>100,
         ExprClass = ifelse(Expr < 100, "no MYCN", ifelse(MYCNStatus != "MNA", "low MYCN", "MYCN-amplified"))) %>%
  ungroup() %>% 
  mutate(CellType = forcats::fct_reorder(CellType, Expr, mean)) %>% 
  mutate(ExprClass = factor(ExprClass, levels=c("no MYCN", "low MYCN", "MYCN-amplified")))


# H3K27ac Fold Change for non-MNA by cell line
as_tibble(h3k27ac_FC) %>%
  mutate(CellType = CellType.x) %>% dplyr::select(-CellType.y, -CellType.x) %>% 
  inner_join(expression_overview, by="CellType") %>%  
  mutate(CellType = forcats::fct_reorder(CellType, Expr, mean)) %>% 
  mutate(bin = paste0(seqnames, "_", start, "_", end)) %>%
  group_by(CellType, ExprClass, bin, seqnames, start, end) %>%
  summarise(FC = mean(FC)) %>%
  ungroup() %>%
  filter(ExprClass != "MYCN-amplified") %>%
  ggplot(aes(x=start, y=FC, color=ExprClass)) +
  geom_col() + 
  facet_grid(CellType ~ .) +
  annotate(xmin=15307032, xmax=15701454, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # NBAS
  annotate(xmin=15731302, xmax=15771235, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # DDX1
  annotate(xmin=16060521, xmax=16076139, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # MYCNUT
  annotate(xmin=16080683, xmax=16087129, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # MYCN
  annotate(xmin=16190549, xmax=16225923, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # GACAT3
  annotate(xmin=16730727, xmax=16847599, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # FAM49A
  theme_kons2() + 
  scale_color_manual(values=c("no MYCN"="steelblue", "MYCN-amplified"="black", "low MYCN"="firebrick3")) +
  theme(strip.text.y = element_text(angle = 0, size=8)) + 
  xlim(gene_of_interest_start-500000,gene_of_interest_end + 500000) + 
  scale_y_continuous(breaks = c(0,30), limits=c(0,30), expand=c(0,0)) +
  theme(strip.text = element_text(size = 8, face = "plain")) +
  theme(panel.spacing.y = unit(6, units="pt")) +
  xlab("") + ylab("H3K7ac [Fold Change over Input]") +
  guides(color=F) +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/BoevaData_H3K27acFC_nonMNA.pdf", 
         height = 3, width = 3,
         useDingbats=F)

# H3K27ac Fold Change for all cell lines
as_tibble(h3k27ac) %>%
  #mutate(CellType = CellType.x) %>% dplyr::select(-CellType.y, -CellType.x) %>% 
  inner_join(expression_overview, by="CellType") %>%  
  mutate(CellType = forcats::fct_reorder(CellType, Expr, mean)) %>% 
  mutate(bin = paste0(seqnames, "_", start, "_", end)) %>%
  group_by(CellType, ExprClass, bin, seqnames, start, end) %>%
  summarise(h3k27ac = mean(h3k27ac)) %>%
  ungroup() %>%
  ggplot(aes(x=start, y=h3k27ac, color=ExprClass)) +
  geom_col() + 
  facet_grid(CellType ~ ., scales = "free_y") +
  annotate(xmin=15307032, xmax=15701454, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # NBAS
  annotate(xmin=15731302, xmax=15771235, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # DDX1
  annotate(xmin=16060521, xmax=16076139, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # MYCNUT
  annotate(xmin=16080683, xmax=16087129, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # MYCN
  annotate(xmin=16190549, xmax=16225923, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # GACAT3
  annotate(xmin=16730727, xmax=16847599, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # FAM49A
  theme_kons2() + 
  scale_color_manual(values=c("no MYCN"="steelblue", "MYCN-amplified"="black", "low MYCN"="firebrick3")) +
  theme(strip.text.y = element_text(angle = 0, size=8)) + 
  xlim(gene_of_interest_start-1000000,gene_of_interest_end + 1000000) + 
  #scale_y_continuous(breaks = c(0,30), limits=c(0,30), expand=c(0,0)) +
  theme(strip.text = element_text(size = 8, face = "plain")) +
  theme(panel.spacing.y = unit(6, units="pt")) +
  xlab("") + ylab("H3K7ac") +
  guides(color=F) +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/BoevaData_H3K27ac_allNBCellLines.pdf", 
         height = 10, width = 3,
         useDingbats=F)


# H3K27ac Fold Change for non-MNA, mean over cell lines
mFC.fig = as_tibble(h3k27ac_FC) %>%
  mutate(CellType = CellType.x) %>% dplyr::select(-CellType.y, -CellType.x) %>% 
  inner_join(expression_overview, by="CellType") %>%  
  mutate(CellType = forcats::fct_reorder(CellType, Expr, mean)) %>% 
  mutate(bin = paste0(seqnames, "_", start, "_", end)) %>%
  group_by(ExprClass, bin, seqnames, start, end) %>%
  summarise(FC = mean(FC)) %>%
  ungroup() %>%
  filter(ExprClass != "MYCN-amplified") %>%
  ggplot(aes(x=start, y=FC, color=ExprClass)) +
  geom_col() + 
  facet_grid(ExprClass ~ .) +
  # annotate(xmin=15307032, xmax=15701454, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # NBAS
  # annotate(xmin=15731302, xmax=15771235, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # DDX1
  # annotate(xmin=16060521, xmax=16076139, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # MYCNUT
  # annotate(xmin=16080683, xmax=16087129, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # MYCN
  # annotate(xmin=16190549, xmax=16225923, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # GACAT3
  # annotate(xmin=16730727, xmax=16847599, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # FAM49A
  # annotate(xmin=cRE$start, xmax=cRE$end, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="red", alpha=0.2) + # cRE
  theme_kons2() + 
  scale_color_manual(values=c("no MYCN"="steelblue", "MYCN-amplified"="black", "low MYCN"="firebrick3")) +
  theme(strip.text.y = element_text(angle = 0, size=8)) + 
  #      axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  xlim(gene_of_interest_start-500000,gene_of_interest_end + 500000) + 
  scale_y_continuous(breaks = c(0,15), limits=c(0,15), expand=c(0,0)) +
  theme(strip.text = element_text(size = 8, face = "plain")) +
  theme(panel.spacing.y = unit(6, units="pt")) +
  xlab("") + ylab("H3K7ac [Fold Change over Input]") +
  guides(color=F) #+
  # ggsave("/Volumes/Elements/MYCNAmplicon/Results/BoevaData_H3K27acFC_nonMNA_MeanFC.pdf", 
  #        height = 1, width = 3,
  #        useDingbats=F)

# H3K27ac Fold Change for non-MNA, DeltaFC
dmFC.fig = as_tibble(h3k27ac_FC) %>%
  mutate(CellType = CellType.x) %>% dplyr::select(-CellType.y, -CellType.x) %>% 
  inner_join(expression_overview, by="CellType") %>%  
  mutate(CellType = forcats::fct_reorder(CellType, Expr, mean)) %>% 
  mutate(bin = paste0(seqnames, "_", start, "_", end)) %>%
  group_by(ExprClass, bin, seqnames, start, end) %>%
  summarise(FC = mean(FC)) %>%
  ungroup() %>%
  filter(ExprClass != "MYCN-amplified") %>%
  tidyr::spread(ExprClass, FC) %>% 
  mutate(DeltaFC = `low MYCN` - `no MYCN`) %>% 
  ggplot(aes(x=start, y=DeltaFC)) +
  geom_col(color="black", fill="black") + 
  # annotate(xmin=15307032, xmax=15701454, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # NBAS
  # annotate(xmin=15731302, xmax=15771235, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # DDX1
  # annotate(xmin=16060521, xmax=16076139, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # MYCNUT
  # annotate(xmin=16080683, xmax=16087129, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # MYCN
  # annotate(xmin=16190549, xmax=16225923, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # GACAT3
  # annotate(xmin=16730727, xmax=16847599, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # FAM49A
  # annotate(xmin=cRE$start, xmax=cRE$end, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="red", alpha=0.2) + # cRE
  theme_kons2() + 
  xlim(gene_of_interest_start-500000,gene_of_interest_end + 500000) + 
  scale_y_continuous(breaks = c(-5,0,15), limits=c(-5,15), expand=c(0,0)) +
  theme(strip.text = element_text(size = 8, face = "plain")) +
  theme(panel.spacing.y = unit(6, units="pt")) +  xlab("") + ylab("H3K7ac Delta FC") 


CLBGA_GATA3_bw = read_bigwig("/Volumes/Elements/nb-cl-chipseq-results/bam/Boeva_CLB-GA_GATA3.trimmed.bwa_hg19.rmdup.filtered.bw", overlap_ranges = roi)
seqlevels(CLBGA_GATA3_bw , pruning.mode="coarse") = standardchrs
CLBGA_GATA3_bw = binnedAverage(bins, coverage(CLBGA_GATA3_bw, weight=CLBGA_GATA3_bw$score), "score")
CLBGA_GATA3_bw = CLBGA_GATA3_bw %>% as_tibble() %>% mutate(Sample = "CLBGA_GATA3")

CLBGA_PHOX2B_bw = read_bigwig("/Volumes/Elements/nb-cl-chipseq-results/bam/Boeva_CLB-GA_PHOX2B.trimmed.bwa_hg19.rmdup.filtered.bw", overlap_ranges = roi)
CLBGA_PHOX2B_bw = binnedAverage(bins, coverage(CLBGA_PHOX2B_bw, weight=CLBGA_PHOX2B_bw$score), "score")
CLBGA_PHOX2B_bw = CLBGA_PHOX2B_bw %>% as_tibble() %>% mutate(Sample = "CLBGA_PHOX2B")
CLBGA_PHOX2B_bw = CLBGA_PHOX2B_bw %>% as_tibble() %>% mutate(Sample = "CLBGA_PHOX2B")

CLBGA_HAND2_bw = read_bigwig("/Volumes/Elements/nb-cl-chipseq-results/bam/Boeva_CLB-GA_HAND2.trimmed.bwa_hg19.rmdup.filtered.bw", overlap_ranges = roi)
CLBGA_HAND2_bw = binnedAverage(bins, coverage(CLBGA_HAND2_bw, weight=CLBGA_HAND2_bw$score), "score")
CLBGA_HAND2_bw = CLBGA_HAND2_bw %>% as_tibble() %>% mutate(Sample = "CLBGA_HAND2")

CLBGA_TF = bind_rows(CLBGA_HAND2_bw, CLBGA_GATA3_bw, CLBGA_PHOX2B_bw)


CLBGA_GATA3.fig = as_tibble(CLBGA_GATA3_bw) %>%
  ggplot(aes(x=start, y=score)) +
  geom_line(color="#238b45") + 
  facet_grid(Sample~.)+
  # annotate(xmin=15307032, xmax=15701454, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # NBAS
  # annotate(xmin=15731302, xmax=15771235, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # DDX1
  # annotate(xmin=16060521, xmax=16076139, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # MYCNUT
  # annotate(xmin=16080683, xmax=16087129, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # MYCN
  # annotate(xmin=16190549, xmax=16225923, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # GACAT3
  # annotate(xmin=16730727, xmax=16847599, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # FAM49A
  # annotate(xmin=cRE$start, xmax=cRE$end, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="red", alpha=0.2) + # cRE
  theme_kons2() + 
  xlim(gene_of_interest_start-500000,gene_of_interest_end + 500000) + 
  scale_y_continuous(breaks = c(0,2), limits=c(0,2), expand=c(0,0)) +
  theme(strip.text = element_text(size = 8, face = "plain")) +
  theme(panel.spacing.y = unit(6, units="pt")) +  xlab("") + 
  ylab("CPM") 

CLBGA_HAND2.fig = as_tibble(CLBGA_HAND2_bw) %>%
  ggplot(aes(x=start, y=score)) +
  geom_line(color="#238b45") + 
  facet_grid(Sample~.)+
  # annotate(xmin=15307032, xmax=15701454, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # NBAS
  # annotate(xmin=15731302, xmax=15771235, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # DDX1
  # annotate(xmin=16060521, xmax=16076139, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # MYCNUT
  # annotate(xmin=16080683, xmax=16087129, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # MYCN
  # annotate(xmin=16190549, xmax=16225923, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # GACAT3
  # annotate(xmin=16730727, xmax=16847599, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # FAM49A
  # annotate(xmin=cRE$start, xmax=cRE$end, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="red", alpha=0.2) + # cRE
  theme_kons2() + 
  xlim(gene_of_interest_start-500000,gene_of_interest_end + 500000) + 
  scale_y_continuous(breaks = c(0,10), limits=c(0,10), expand=c(0,0)) +
  theme(strip.text = element_text(size = 8, face = "plain")) +
  theme(panel.spacing.y = unit(6, units="pt")) +  xlab("") + 
  ylab("CPM") 

CLBGA_PHOX2B.fig = as_tibble(CLBGA_PHOX2B_bw) %>%
  ggplot(aes(x=start, y=score)) +
  geom_col(color="#238b45") + 
  facet_grid(Sample~.)+
  # annotate(xmin=15307032, xmax=15701454, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # NBAS
  # annotate(xmin=15731302, xmax=15771235, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # DDX1
  # annotate(xmin=16060521, xmax=16076139, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # MYCNUT
  # annotate(xmin=16080683, xmax=16087129, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # MYCN
  # annotate(xmin=16190549, xmax=16225923, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # GACAT3
  # annotate(xmin=16730727, xmax=16847599, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # FAM49A
  # annotate(xmin=cRE$start, xmax=cRE$end, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="red", alpha=0.2) + # cRE
  theme_kons2() + 
  xlim(gene_of_interest_start-500000,gene_of_interest_end + 500000) + 
  scale_y_continuous(breaks = c(0,10), limits=c(0,10), expand=c(0,0)) +
  theme(strip.text = element_text(size = 8, face = "plain")) +
  theme(panel.spacing.y = unit(6, units="pt")) +  xlab("") + 
  ylab("CPM") 

CLBGA_TF.fig = as_tibble(CLBGA_TF) %>%
  ggplot(aes(x=start, y=score)) +
  geom_col(color="#238b45") + 
  facet_grid(Sample~.)+
  # annotate(xmin=15307032, xmax=15701454, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # NBAS
  # annotate(xmin=15731302, xmax=15771235, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # DDX1
  # annotate(xmin=16060521, xmax=16076139, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # MYCNUT
  # annotate(xmin=16080683, xmax=16087129, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # MYCN
  # annotate(xmin=16190549, xmax=16225923, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # GACAT3
  # annotate(xmin=16730727, xmax=16847599, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # FAM49A
  # annotate(xmin=genes_df$start, xmax=genes_df$end, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="black", alpha=0.2) + # genes
  # annotate(xmin=cRE$start, xmax=cRE$end, ymin=-Inf, ymax=Inf, geom="rect", color=NA, fill="red", alpha=0.2) + # cRE
  theme_kons2() + 
  xlim(gene_of_interest_start-500000,gene_of_interest_end + 500000) + 
  scale_y_continuous(breaks = c(0,10), limits=c(0,10), expand=c(0,0)) +
  theme(strip.text = element_text(size = 8, face = "plain")) +
  theme(panel.spacing.y = unit(6, units="pt")) +  xlab("") + 
  ylab("CPM") 

genes_df = 
  data.frame(
    "start" = c(15307032, 15731302, 16060521, 16080683, 16190549, 16730727),
    "end" = c(15701454, 15771235, 16076139, 16087129, 16225923, 16847599),
    "name" = c("NBAS", "DDX1", "MYCNUT", "MYCN", "GACAT3", "FAM49A"),
    "class" = rep("gene", 6)
)

genes.fig = genes_df %>%
  ggplot(aes(x=start, y=1)) +
  geom_rect(xmin=genes_df$start, xmax=genes_df$end, ymin=-Inf, ymax=Inf, color=NA, fill="black") + 
  theme_kons2() + 
  theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank()) +
  xlim(gene_of_interest_start-500000,gene_of_interest_end + 500000)

cRE.fig = cRE %>%
  ggplot(aes(x=start, y=maxDFC)) +
  geom_rect(xmin=cRE$start, xmax=cRE$end, ymin=-Inf, ymax=Inf, color="firebrick3", fill="firebrick3") + 
  theme_kons2() + 
  theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank()) +
  xlim(gene_of_interest_start-500000,gene_of_interest_end + 500000)

cRE$class = "cRE"
cRE_genes_plot_df = 
  bind_rows(cRE, genes_df) %>% 
  mutate(dummy=1) %>% 
  dplyr::select(start, end, class, dummy) %>%
  mutate(this_color = ifelse(class == "gene", NA, "firebrick3"),
         this_fill = ifelse(class == "gene", "black", "firebrick3"))

cRE_genes.fig = 
  cRE_genes_plot_df %>%
  ggplot(aes(x=start, y=dummy)) +
  geom_rect(xmin=cRE_genes_plot_df$start, xmax=cRE_genes_plot_df$end, ymin=-Inf, ymax=Inf, color=cRE_genes_plot_df$this_color, fill=cRE_genes_plot_df$this_fill) + 
  theme_kons2() +
  theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank()) +
  xlim(gene_of_interest_start-500000,gene_of_interest_end + 500000)

# ggsave("/Volumes/Elements/MYCNAmplicon/Results/BoevaData_H3K27acFC_nonMNA_MeanAndDeltaFC.pdf",
#        egg::ggarrange(mFC.fig + 
#                  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
#                dmFC.fig +
#                  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
#                CLBGA_TF.fig,
#                nrow = 3),
#        height = 3, width=5, onefile = FALSE, useDingbats = F)

ggsave("/Volumes/Elements/MYCNAmplicon/Results/BoevaData_H3K27acFC_nonMNA_DeltaFCAndTF_SingleTFPlots.pdf",
       egg::ggarrange(
         genes.fig +
           theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
         dmFC.fig +
           theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
         cRE.fig +
           theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
         CLBGA_PHOX2B.fig+
           theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
         CLBGA_GATA3.fig+
           theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
         CLBGA_HAND2.fig,
         nrow = 6, 
         heights = c(0.25,2,0.25,1,1,1)),
       height = 3, width=3, onefile = FALSE, useDingbats = F)

# ------------------------------------------------------------------------------
# Generate genome-wide mean H3K27ac tracks for different classes of genes
# ------------------------------------------------------------------------------

metadata_chipseq = 
  metadata_chipseq %>% filter(grepl("Boeva", H3K27ac_bw), !grepl("hNCC", H3K27ac_bw))

standardchrs = standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)[1:23]
binsize = 1000 
bins = tileGenome(seqlengths = seqlengths(BSgenome.Hsapiens.UCSC.hg19),tilewidth = binsize, cut.last.tile.in.chrom = T)
seqlevels(bins, pruning.mode="coarse") = standardchrs

controls = list()
h3k27ac = list()
h3k27ac_FC = list()
h3k27ac_peaks = list()
h3k27ac_peaks_unbinned = list()
for (i in 1:nrow(metadata_chipseq)){
  print(i)
  
  controls[[i]] = read_bigwig(metadata_chipseq[[i,"Input_bw"]])#, overlap_ranges = roi)
  seqlevels(controls[[i]] , pruning.mode="coarse") = standardchrs
  controls[[i]] = binnedAverage(bins, coverage(controls[[i]], weight=controls[[i]]$score), "input")
  controls[[i]]$CellType = metadata_chipseq[[i,"CellType"]]
  
  h3k27ac[[i]] = read_bigwig(metadata_chipseq[[i,"H3K27ac_bw"]])#, overlap_ranges = roi)
  seqlevels(h3k27ac[[i]] , pruning.mode="coarse") = standardchrs
  h3k27ac[[i]] = binnedAverage(bins, coverage(h3k27ac[[i]], weight=h3k27ac[[i]]$score), "h3k27ac")
  h3k27ac[[i]]$CellType = metadata_chipseq[[i,"CellType"]]
  
  h3k27ac_FC[[i]] = 
    join_overlap_inner(controls[[i]], h3k27ac[[i]], minoverlap=binsize) %>% 
    mutate(FC = (h3k27ac+0.01)/(input+0.01))
  
  h3k27ac_peaks[[i]] = read_narrowpeaks(metadata_chipseq[[i,"H3K27ac_peaks"]])#), overlap_ranges = roi)
  seqlevels(h3k27ac_peaks[[i]] , pruning.mode="coarse") = standardchrs
  h3k27ac_peaks[[i]] = binnedAverage(bins, coverage(h3k27ac_peaks[[i]]), "hasPeaks") %>% mutate(hasPeaks = (hasPeaks>0))
  h3k27ac_peaks[[i]]$CellType = metadata_chipseq[[i,"CellType"]]
}
controls = do.call(c,controls)
h3k27ac = do.call(c,h3k27ac)
h3k27ac_FC = do.call(c,h3k27ac_FC)
h3k27ac_peaks = do.call(c,h3k27ac_peaks)

#save.image("/Volumes/Elements/MYCNAmplicon/Results/BoevaH3K27acData_AllChr.Rdata")
#load("/Volumes/Elements/MYCNAmplicon/Results/BoevaH3K27acData_AllChr.Rdata")

# Compute genome-wide mean fold change for non-MYCN-expressing cell lines (7 cell lines)
MeanFC_nonMYCNExpr = as_tibble(h3k27ac_FC) %>%
  filter(CellType.x %in% c("GICAN", "SHEP", "SKNAS", "GIMEN", "SKNSH", "NB69", "SJNB12")) %>%
  group_by(seqnames, start, end) %>%
  summarise(score = mean(FC, na.rm=T)) 
MeanFC_nonMYCNExpr = makeGRangesFromDataFrame(MeanFC_nonMYCNExpr, keep.extra.columns = T)
seqlevels(MeanFC_nonMYCNExpr, pruning.mode="coarse") = seqlevels(BSgenome.Hsapiens.UCSC.hg19)
seqlengths(MeanFC_nonMYCNExpr) = seqlengths(BSgenome.Hsapiens.UCSC.hg19)
MeanFC_nonMYCNExpr %>%
  write_bigwig("/Volumes/Elements/MYCNAmplicon/Results/Boeva_nonMYCNExpr_H3K27ac_MeanFC.bw")

# Compute genome-wide mean fold change for MYCN-expressing, MYCN-non-amplified cell lines (5 cell lines)
MeanFC_lowMYCNExpr = as_tibble(h3k27ac_FC) %>%
  filter(CellType.x %in% c("SHSY5Y", "SJNB1", "SKNFI", "CLBGA", "NBEBc1")) %>% 
  group_by(seqnames, start, end) %>%
  summarise(score = mean(FC, na.rm=T))
MeanFC_lowMYCNExpr = makeGRangesFromDataFrame(MeanFC_lowMYCNExpr, keep.extra.columns = T)
seqlevels(MeanFC_lowMYCNExpr, pruning.mode="coarse") = seqlevels(BSgenome.Hsapiens.UCSC.hg19)
seqlengths(MeanFC_lowMYCNExpr) = seqlengths(BSgenome.Hsapiens.UCSC.hg19)
MeanFC_lowMYCNExpr %>%
  write_bigwig("/Volumes/Elements/MYCNAmplicon/Results/Boeva_lowMYCNExpr_H3K27ac_MeanFC.bw")

# Compute genome-wide mean fold change for MYCN-amplified cell lines (13 cell lines)
MeanFC_MNA = as_tibble(h3k27ac_FC) %>%
  filter(CellType.x %in% c("LAN1", "CLBPE", "SKNDZ", "CLBCAR", "CLBMA",
                           "IMR32", "CHP212", "SJNB8", "TR14", "SKNBE2C",
                           "N206", "SJNB6", "CLBBERLud")) %>% 
  group_by(seqnames, start, end) %>%
  dplyr::summarise(score = mean(FC, na.rm=T))
MeanFC_MNA = makeGRangesFromDataFrame(MeanFC_MNA, keep.extra.columns = T)
seqlevels(MeanFC_MNA, pruning.mode="coarse") = seqlevels(BSgenome.Hsapiens.UCSC.hg19)
seqlengths(MeanFC_MNA) = seqlengths(BSgenome.Hsapiens.UCSC.hg19)
MeanFC_MNA %>%
  write_bigwig("/Volumes/Elements/MYCNAmplicon/Results/Boeva_MNA_H3K27ac_MeanFC.bw")

# Calculate genome-wide Delta FC: lowMYCNExpr - nonMYCNExpr
MeanFC_nonMYCNExpr = MeanFC_nonMYCNExpr %>%
  mutate(nonMYCN_meanFC = score) %>%
  dplyr::select(seqnames, start, end, nonMYCN_meanFC)
MeanFC_lowMYCNExpr = MeanFC_lowMYCNExpr %>%
  mutate(lowMYCN_meanFC = score) %>%
  dplyr::select(seqnames, start, end, lowMYCN_meanFC)
DeltaFC_LowMinusNon = full_join(MeanFC_lowMYCNExpr, MeanFC_nonMYCNExpr) %>%
  mutate(score = lowMYCN_meanFC-nonMYCN_meanFC) %>%
  dplyr::select(seqnames, start, end, score)
DeltaFC_LowMinusNon = makeGRangesFromDataFrame(DeltaFC_LowMinusNon, keep.extra.columns = T)
seqlevels(DeltaFC_LowMinusNon, pruning.mode="coarse") = seqlevels(BSgenome.Hsapiens.UCSC.hg19)
seqlengths(DeltaFC_LowMinusNon) = seqlengths(BSgenome.Hsapiens.UCSC.hg19)
DeltaFC_LowMinusNon %>%
  write_bigwig("/Volumes/Elements/MYCNAmplicon/Results/Boeva_lowMinusNon_H3K27ac_DeltaMeanFC.bw")

# ------------------------------------------------------------------------------
# bla
# ------------------------------------------------------------------------------

metadata_chipseq = 
  metadata_chipseq %>% filter(grepl("Boeva", H3K27ac_bw), !grepl("hNCC", H3K27ac_bw))

standardchrs = standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)[1:23]
binsize = 1000 
bins = tileGenome(seqlengths = seqlengths(BSgenome.Hsapiens.UCSC.hg19),tilewidth = binsize, cut.last.tile.in.chrom = T)
seqlevels(bins, pruning.mode="coarse") = standardchrs

# Read peaks files
h3k27ac_peaks_unbinned = list()
for (i in 1:nrow(metadata_chipseq)){
  print(i)
  h3k27ac_peaks_unbinned[[i]] = read_narrowpeaks(metadata_chipseq[[i,"H3K27ac_peaks"]])
  seqlevels(h3k27ac_peaks_unbinned[[i]] , pruning.mode="coarse") = standardchrs
  h3k27ac_peaks_unbinned[[i]]$CellType = metadata_chipseq[[i,"CellType"]]
}
h3k27ac_peaks_unbinned = do.call(c,h3k27ac_peaks_unbinned)

# Peaks for low MYCN expression
peak_regions_allchr = 
  coverage(h3k27ac_peaks_unbinned %>% 
             filter(CellType %in% mycnexp_non_mna_celltypes)) %>%
  GRanges() %>%
  filter(score >= 3) %>% 
  reduce(min.gapwidth = 2000)
seqlevels(peak_regions_allchr , pruning.mode="coarse") = standardchrs
seqlengths(peak_regions_allchr) = seqlengths(BSgenome.Hsapiens.UCSC.hg19)[standardchrs]
peak_regions_allchr %>%
  write_bed(file = "/Volumes/Elements/MYCNAmplicon/Results/Boeva_lowMYCNExpr_H3K27ac_AggregatePeaks.bed")

# Peaks for non-MYCN expression
peak_regions_allchr_nonexpr = 
  coverage(h3k27ac_peaks_unbinned %>% 
             filter(CellType %in% c("GICAN", "SHEP", "SKNAS", "GIMEN", "SKNSH", "NB69", "SJNB12"))) %>%
  GRanges() %>%
  filter(score >= 4) %>% 
  reduce(min.gapwidth = 2000)
seqlevels(peak_regions_allchr_nonexpr , pruning.mode="coarse") = standardchrs
seqlengths(peak_regions_allchr_nonexpr) = seqlengths(BSgenome.Hsapiens.UCSC.hg19)[standardchrs]

peak_regions_allchr_nonexpr %>%
  write_bed(file = "/Volumes/Elements/MYCNAmplicon/Results/Boeva_noMYCNExpr_H3K27ac_AggregatePeaks.bed")

CLBGA_GATA3_peaks = read_narrowpeaks("/Volumes/Elements/nb-cl-chipseq-results/MACS2/Boeva_CLB-GA_GATA3/Boeva_CLB-GA_GATA3_MACS2_nocontrol_peaks.narrowPeak")
seqlevels(CLBGA_GATA3_peaks , pruning.mode="coarse") = standardchrs
seqlengths(CLBGA_GATA3_peaks) = seqlengths(BSgenome.Hsapiens.UCSC.hg19)[standardchrs]
CLBGA_PHOX2B_peaks = read_narrowpeaks("/Volumes/Elements/nb-cl-chipseq-results/MACS2/Boeva_CLB-GA_PHOX2B/Boeva_CLB-GA_PHOX2B_MACS2_nocontrol_peaks.narrowPeak")
seqlevels(CLBGA_PHOX2B_peaks , pruning.mode="coarse") = standardchrs
seqlengths(CLBGA_PHOX2B_peaks) = seqlengths(BSgenome.Hsapiens.UCSC.hg19)[standardchrs]
CLBGA_HAND2_peaks = read_narrowpeaks("/Volumes/Elements/nb-cl-chipseq-results/MACS2/Boeva_CLB-GA_HAND2/Boeva_CLB-GA_HAND2_MACS2_nocontrol_peaks.narrowPeak")
seqlevels(CLBGA_HAND2_peaks , pruning.mode="coarse") = standardchrs
seqlengths(CLBGA_HAND2_peaks) = seqlengths(BSgenome.Hsapiens.UCSC.hg19)[standardchrs]
CLBGA_CRC_peaks = GRanges(coverage(c(CLBGA_GATA3_peaks, CLBGA_PHOX2B_peaks, CLBGA_HAND2_peaks))) %>% filter(score > 0)
seqlevels(CLBGA_CRC_peaks , pruning.mode="coarse") = standardchrs

peak_regions_allchr %>%
  filter_by_overlaps(CLBGA_CRC_peaks) %>% 
  write_bed(file = "/Volumes/Elements/MYCNAmplicon/Results/Boeva_lowMYCNExpr_H3K27ac_AggregatePeaks_CRCFactorPositiveOnly.bed")

# ------------------------------------------------------------------------------
# how large are the enhancers and how far are they away from MYCN
# ------------------------------------------------------------------------------

enhancers = read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_nMNA_MYCN_Enhancers.bed")

# size
enhancer_sizes = width(enhancers)
names(enhancer_sizes) = enhancers$name
print(enhancer_sizes)

# distance to the MYCN TSS
mycn_tss = 16080683
enhancer_distance_to_mycn = pmin(abs(start(enhancers)-mycn_tss), abs(end(enhancers)-mycn_tss))
names(enhancer_distance_to_mycn) = enhancers$name
print(enhancer_distance_to_mycn)


# ------------------------------------------------------------------------------
# generate aggregate SE
# ------------------------------------------------------------------------------

noMYCNSamples = c("GICAN", "SHEP", "SKNAS", "GIMEN", "SKNSH", "NB69", "SJNB12")
lowMYCNSamples = c("SHSY5Y", "SJNB1", "SKNFI", "CLBGA", "NBEBc1")
MNASamples = c("LAN1", "CLBPE", "SKNDZ", "CLBCAR", "CLBMA",
               "IMR32", "CHP212", "SJNB8", "TR14", "SKNBE2C",
               "N206", "SJNB6", "CLBBERLud")

metadata_chipseq = 
  metadata_chipseq %>% filter(grepl("Boeva", H3K27ac_bw), !grepl("hNCC", H3K27ac_bw))
metadata_chipseq$Sample = metadata_chipseq$H3K27ac_bw %>%
  gsub("/Volumes/Elements/nb-cl-chipseq-results/bam/", "", .) %>%
  gsub(".trimmed.bwa_hg19.rmdup.filtered.bw", "", .)
metadata_chipseq = metadata_chipseq %>% dplyr::select(CellType, Sample)
metadata_chipseq$lily_fname = paste0("/Volumes/Elements/nb-cl-chipseq-results/LILY/",
                                     metadata_chipseq$Sample,
                                     "_hmcan.scores.bed")
metadata_chipseq$ExprClass = 
  ifelse(metadata_chipseq$CellType %in% noMYCNSamples,
         "noMYCN",
         ifelse(
           metadata_chipseq$CellType %in% lowMYCNSamples,
           "lowMYCN", 
           ifelse(metadata_chipseq$CellType %in% MNASamples, 
           "MNA", 
           NA
         )))

lily = lapply(1:nrow(metadata_chipseq), function(i) read_bed(metadata_chipseq[[i,"lily_fname"]]) %>% mutate(Sample = metadata_chipseq[[i, "CellType"]],
                                                                                                            ExprClass = metadata_chipseq[[i, "ExprClass"]]))
lily = do.call(c, lily)
lily = lily %>% filter(name == "SE")

noMYCNSE = lily %>% 
  filter(ExprClass == "noMYCN") %>%
  reduce()
noMYCNSE %>%
  write_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_noMYCN_SE.bed")

lowMYCNSE = lily %>% 
  filter(ExprClass == "lowMYCN") %>%
  reduce()
lowMYCNSE %>%
  write_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_lowMYCN_SE.bed")

lily %>% 
  filter(ExprClass == "lowMYCN") %>%
  coverage() %>%
  GRanges() %>%
  filter(score >= 2) %>%
  reduce() %>%
  as_tibble() %>% 
  View

nonMNA_SE = lily %>%
  filter(ExprClass %in% c("lowMYCN", "noMYCN")) %>%
  reduce()
nonMNA_SE %>%
  write_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_noMNA_SE.bed")

CLBGA_CRC_peaks = GRanges(coverage(c(CLBGA_GATA3_peaks, CLBGA_PHOX2B_peaks, CLBGA_HAND2_peaks))) %>% filter(score > 0)
seqlevels(CLBGA_CRC_peaks , pruning.mode="coarse") = standardchrs
lowMYCN_CRCdriven_SE = 
  lowMYCNSE %>% filter_by_overlaps(CLBGA_CRC_peaks)
lowMYCN_CRCdriven_SE %>%
  write_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_lowMYCN_CRCdriven_SE.bed")
