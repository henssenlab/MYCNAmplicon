library(dplyr)
library(plyranges)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg19)
library(biomaRt)
source("/Volumes/Elements/MYCNAmplicon/Code/CustomThemes.R")
source("/Volumes/Elements/MYCNAmplicon/Code/hasOverlap.R")

gene_of_interest = "MYCN"
gene_of_interest_ensembl = "ENSG00000134323"
gene_of_interest_chr = "chr2"
gene_of_interest_start = 16080683
gene_of_interest_end = 16087129
window_of_interest = 1200000
roi = GRanges(
  seqnames = gene_of_interest_chr,
  ranges = IRanges(start = gene_of_interest_start - window_of_interest,
                   end = gene_of_interest_end + window_of_interest),
  strand = "*"
)
roi = GRanges(
  seqnames = gene_of_interest_chr,
  ranges = IRanges(start = 14600000,
                   end = 17500000),
  strand = "*"
)
standardchrs = standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)[1:23]
binsize = 1000 
bins = tileGenome(seqlengths = seqlengths(BSgenome.Hsapiens.UCSC.hg19),tilewidth = binsize, cut.last.tile.in.chrom = T) %>%
  filter_by_overlaps(roi, minoverlap = binsize/10)
seqlevels(bins, pruning.mode="coarse") = standardchrs

make_figure = function(fname, roi, plot.y.min=0, plot.y.max=100, this_color = "black", ftype = "bigwig"){
  
  if (ftype == "bigwig"){
    this_track = read_bigwig(fname, overlap_ranges = roi) 
  } else if (ftype == "bedgraph"){
    this_track = read_bed_graph(fname, overlap_ranges = roi) 
  } else {
    stop("Invalid ftype.")
  }

  standardchrs = standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)[1:23]
  binsize = 1000 
  bins = tileGenome(seqlengths = seqlengths(BSgenome.Hsapiens.UCSC.hg19),tilewidth = binsize, cut.last.tile.in.chrom = T) %>%
    filter_by_overlaps(roi, minoverlap = binsize/10)
  seqlevels(bins, pruning.mode="coarse") = standardchrs
  
  seqlevels(this_track , pruning.mode="coarse") = standardchrs
  this_track = binnedAverage(bins, coverage(this_track, weight=this_track$score), "score")
  fig = 
    this_track %>%
    as_tibble() %>% 
    ggplot(aes(x=start, y=score)) +
    geom_area(fill=this_color, color=this_color) + 
    theme_kons2() +
    ylab("") +
    scale_y_continuous(breaks = c(plot.y.min,plot.y.max), limits=c(plot.y.min,plot.y.max), expand=c(0,0)) 
  return(fig)
}
round_up_to = function(s, stepsize=50) ceiling(s/stepsize)*stepsize

make_copynumber_figure = function(fname, roi){
  #fname = "/Volumes/Elements/nb-cl-wgs/nb-cl-wgs-freec/Kelly/Kelly_IlluminaWGS.hg19.rmdup.fixRG.bam_CNVs.fillgaps.txt"
  #fname = "/Volumes/Elements/nb-cl-wgs/nb-cl-wgs-freec/NGP/NGP_IlluminaWGS.hg19.rmdup.fixRG.bam_CNVs.fillgaps.txt"
  copynumber_data = read.table(fname) %>%
    dplyr::select(V1, V2, V3, V6, V7)
  colnames(copynumber_data) = c("seqnames", "start", "end", "copynumber", "type")
  copynumber_data = 
    copynumber_data %>%
    as_tibble() %>% 
    filter(seqnames == "chr2") %>%
    filter(hasOverlap(start(roi), end(roi), start, end)) %>%
    mutate(start = pmax(start(roi), start),
           end = pmin(end(roi), end))
  fig = 
    copynumber_data %>%
    ggplot() +
    geom_segment(aes(x=start, xend=end, y=copynumber, yend=copynumber), size=0.25) + 
    theme_kons2() +
    xlim(start(roi), end(roi)) +
    ylab("") +
    scale_y_continuous(breaks = c(0,round_up_to(max(copynumber_data$copynumber))), limits=c(0,round_up_to(max(copynumber_data$copynumber))), expand=c(0,0)) 
  return(fig)
}


genes_df = 
  data.frame(
    "start" = c(14772810, 15307032, 15731302, 16060521, 16080683, 16190549, 16730727, 17691851, 17720393, 17845079),
    "end" = c(14790933, 15701454, 15771235, 16076139, 16087129, 16225923, 16847599, 17699706, 17838285, 17981509),
    "name" = c("FAM84A", "NBAS", "DDX1", "MYCNUT", "MYCN", "GACAT3", "FAM49A", "RAD51AP2", "VSNL1", "SMC6"),
    "class" = rep("gene", 10)
  )
genes.fig = genes_df %>%
  ggplot(aes(x=start, y=1)) +
  geom_rect(xmin=genes_df$start, xmax=genes_df$end, ymin=-Inf, ymax=Inf, color=NA, fill="black") + 
  theme_kons2() + 
  theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank()) +
  xlim(start(roi),end(roi))
cRE = read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_nMNA_MYCN_Enhancers.bed")
cRE = as_tibble(cRE)
cRE.fig = cRE %>%
  ggplot(aes(x=start, y=1)) +
  geom_rect(xmin=cRE$start, xmax=cRE$end, ymin=-Inf, ymax=Inf, color="firebrick3", fill="firebrick3") + 
  theme_kons2() + 
  theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank()) +
  xlim(start(roi),end(roi))
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
  xlim(start(roi),end(roi))

v4C_color = RColorBrewer::brewer.pal(4, "Dark2")[1]
ATAC_color = RColorBrewer::brewer.pal(4, "Dark2")[2]
H3K27ac_color = RColorBrewer::brewer.pal(4, "Dark2")[3]
H3K4me1_color = RColorBrewer::brewer.pal(4, "Dark2")[4]
CTCF_color = RColorBrewer::brewer.pal(6, "Dark2")[5]
RAD21_color = RColorBrewer::brewer.pal(6, "Dark2")[6]
meanH3K27ac_color = "gold2" #RColorBrewer::brewer.pal(8, "Dark2")[8]
CRC_color = RColorBrewer::brewer.pal(8, "Dark2")[8]



ggsave(
  "/Volumes/Elements/MYCNAmplicon/Results/Fig2c.pdf",
  egg::ggarrange(
  genes.fig +
    theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
  ggplot() + theme_void(),
  cRE.fig +
    theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
  ggplot() + theme_void(),
  make_copynumber_figure("/Volumes/Elements/nb-cl-wgs/nb-cl-wgs-freec/Kelly/Kelly_IlluminaWGS.hg19.rmdup.fixRG.bam_CNVs.fillgaps.txt", roi=roi) +
    theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
  make_figure("/Volumes/Elements/nb-cl-atacseq-results/bam/Zeid_KELLY_ATAC.trimmed.hg19.rmdup.filterednormed.bw", roi=roi, plot.y.min=0, plot.y.max=35, this_color = ATAC_color) +
    theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
  make_figure("/Volumes/Elements/nb-cl-chipseq-results/bam/MundlosLab_Kelly_H3K27ac.trimmed.bwa_hg19.rmdup.filtered.bw", roi=roi, plot.y.min=0, plot.y.max=35, this_color = H3K27ac_color) +
    theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
  make_figure("/Volumes/Elements/nb-cl-chipseq-results/bam/MundlosLab_Kelly_H3K4me1.trimmed.bwa_hg19.rmdup.filtered.bw", roi=roi, plot.y.min=0, plot.y.max=25, this_color = H3K4me1_color) +
    theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
  make_figure("/Volumes/Elements/MYCNAmplicon/Data/4C/4C_KELLY_L14983.bwa.bam.bedgraph", roi=roi, plot.y.min=0, plot.y.max=18000, this_color=v4C_color, ftype = "bedgraph") +
    theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
  ggplot() + theme_void(),
  make_copynumber_figure("/Volumes/Elements/nb-cl-wgs/nb-cl-wgs-freec/NGP/NGP_IlluminaWGS.hg19.rmdup.fixRG.bam_CNVs.fillgaps.txt", roi=roi) +
    theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
  make_figure("/Volumes/Elements/nb-cl-atacseq-results/bam/Zeid_NGP_ATAC.trimmed.hg19.rmdup.filterednormed.bw", roi=roi, plot.y.min=0, plot.y.max=50, this_color = ATAC_color) +
    theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
  make_figure("/Volumes/Elements/nb-cl-chipseq-results/bam/MundlosLab_NGP_H3K27ac.trimmed.bwa_hg19.rmdup.filtered.bw", roi=roi, plot.y.min=0, plot.y.max=25, this_color = H3K27ac_color) +
    theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
  make_figure("/Volumes/Elements/nb-cl-chipseq-results/bam/MundlosLab_NGP_H3K4me1.trimmed.bwa_hg19.rmdup.filtered.bw", roi=roi, plot.y.min=0, plot.y.max=75, this_color = H3K4me1_color) +
    theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
  make_figure("/Volumes/Elements/MYCNAmplicon/Data/4C/4C_NGP_L14984.bwa.bam.bedgraph", roi=roi, plot.y.min=0, plot.y.max=18000, this_color=v4C_color, ftype = "bedgraph"),
  nrow=15,
  heights = c(0.5,0.5,0.5,1,1,1,1,1,1,1,1,1,1,1,1)),
  height=3.5, width=3, onefile=F, useDingbats=F)
# old was: height=4, width=6
# original 4C height: 100000

# ------------------------------------------------------------------------------
# Plot IMR-5/75 and CHP-212
# ------------------------------------------------------------------------------

ggsave(
  "/Volumes/Elements/MYCNAmplicon/Results/IMR5_CHP_EpigeneticsFig.pdf",
  egg::ggarrange(
    genes.fig +
      theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
    ggplot() + theme_void(),
    cRE.fig +
      theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
    ggplot() + theme_void(),
    make_copynumber_figure("/Volumes/Elements/nb-cl-wgs/nb-cl-wgs-freec/IMR-5-75/IMR-5-75_IlluminaWGS.hg19.rmdup.fixRG.bam_CNVs.fillgaps.txt", roi=roi) +
      theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
  make_figure("/Volumes/Elements/nb-cl-atacseq-results/bam/Mundlos_IMR575_ATAC.trimmed.hg19.rmdup.filterednormed.bw", roi=roi, plot.y.min=0, plot.y.max=50, this_color = ATAC_color) +
      theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
  make_figure("/Volumes/Elements/nb-cl-pe-chipseq-results/bam/Schulte_IMR5_H3K27ac.trimmed.hg19.rmdup.filterednormed.bw", roi=roi, plot.y.min=0, plot.y.max=60, this_color = H3K27ac_color) +
      theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
    make_figure("/Volumes/Elements/nb-cl-chipseq-results/bam/MundlosLab_IMR-5-75_H3K4me1.trimmed.bwa_hg19.rmdup.filtered.bw", roi=roi, plot.y.min=0, plot.y.max=25, this_color = H3K4me1_color) +
      theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
    make_figure("/Volumes/Elements/Virtual4C/IMR5_hg19_canonical_MAPQ30_merged.hic.MYCN4C_Aug28.Normalization_KR_Viewpoint_16075000_16085000_Res_5000.bw", roi=roi, plot.y.min=0, plot.y.max=8, this_color=v4C_color) +
      theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
    ggplot() + theme_void(),
    make_copynumber_figure("/Volumes/Elements/nb-cl-wgs/nb-cl-wgs-freec/CHP-212/CHP-212_IlluminaWGS.hg19.rmdup.fixRG.bam_CNVs.fillgaps.txt", roi=roi) +
      theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
    make_figure("/Volumes/Elements/nb-cl-atacseq-results/bam/Mundlos_CHP212_ATAC.trimmed.hg19.rmdup.filterednormed.bw", roi=roi, plot.y.min=0, plot.y.max=50, this_color = ATAC_color) +
      theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
    make_figure("/Volumes/Elements/nb-cl-chipseq-results/bam/Boeva_CHP212_H3K27ac.trimmed.bwa_hg19.rmdup.filtered.bw", roi=roi, plot.y.min=0, plot.y.max=30, this_color = H3K27ac_color) +
      theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
    make_figure("/Volumes/Elements/nb-cl-chipseq-results/bam/MundlosLab_CHP-212_H3K4me1.trimmed.bwa_hg19.rmdup.filtered.bw", roi=roi, plot.y.min=0, plot.y.max=20, this_color = H3K4me1_color) +
      theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
    make_figure("/Volumes/Elements/Virtual4C/CHP_hg19_canonical_MAPQ30_merged.hic.MYCN4C_Aug28.Normalization_KR_Viewpoint_16075000_16085000_Res_5000.bw", roi=roi, plot.y.min=0, plot.y.max=0.5, this_color=v4C_color),
    nrow=15,
    heights = c(0.5,0.5,0.5,1,1,1,1,1,1,1,1,1,1,1,1)),
  height=3.5, width=3, onefile=F, useDingbats=F)

# ------------------------------------------------------------------------------
# Plot IMR-5/75 fragments
# ------------------------------------------------------------------------------

mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl", 
                   mirror = "useast", 
                   GRCh = 37)

window = GRanges(
  seqnames = rep("chr2", 6),
  ranges = IRanges(start = c(14500000, 17300000, 29500000, 53300000, 66700000, 69100000),
                   end   = c(16200000, 17400000, 30000000, 53700000, 67800000, 69600000))
)

for (r in 1:length(window)){
  roi = window[r]
  
  attributes <- c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "gene_biotype")
  filters <- c("chromosome_name","start","end")
  genes = getBM(attributes=attributes, 
                filters=filters, values=list(chromosome=gsub("chr", "", as.character(seqnames(roi))),
                                                                    start=as.character(start(roi)),
                                                                    end=as.character(end(roi))), 
                mart=mart) %>% 
    as_tibble() %>% 
    filter(hgnc_symbol != "", gene_biotype == "protein_coding") %>%
    mutate(chromosome_name = paste0("chr", chromosome_name)) %>% 
    dplyr::select(hgnc_symbol, chromosome_name, start_position, end_position)
  colnames(genes) = c("gene", "seqnames", "start", "end")
  genes_df = genes
  genes_df = rbind(genes_df, 
                   data.frame(gene="NULLGENE",seqnames="chr2", start=0, end=1))
  #genes_= makeGRangesFromDataFrame(genes, keep.extra.columns = T)
  genes.fig = genes_df %>%
    ggplot(aes(x=start, y=1)) +
    geom_rect(xmin=genes_df$start, xmax=genes_df$end, ymin=-Inf, ymax=Inf, color=NA, fill="black") + 
    theme_kons2() + 
    theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank()) +
    xlim(start(roi),end(roi))

  ggsave(
    paste0("/Volumes/Elements/MYCNAmplicon/Results/IMR5_EpigeneticsFig_Fragments_",
           as.character(seqnames(roi)),"_",as.character(start(roi)),"_", as.character(end(roi)), 
           ".pdf"),
    egg::ggarrange(
      genes.fig +
        theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
      ggplot() + theme_void(),
      make_copynumber_figure("/Volumes/Elements/nb-cl-wgs/nb-cl-wgs-freec/IMR-5-75/IMR-5-75_IlluminaWGS.hg19.rmdup.fixRG.bam_CNVs.fillgaps.txt", roi=roi) +
        theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
      make_figure("/Volumes/Elements/nb-cl-atacseq-results/bam/Mundlos_IMR575_ATAC.trimmed.hg19.rmdup.filterednormed.bw", roi=roi, plot.y.min=0, plot.y.max=50, this_color = ATAC_color) +
        theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
      make_figure("/Volumes/Elements/nb-cl-pe-chipseq-results/bam/Schulte_IMR5_H3K27ac.trimmed.hg19.rmdup.filterednormed.bw", roi=roi, plot.y.min=0, plot.y.max=60, this_color = H3K27ac_color) +
        theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
      make_figure("/Volumes/Elements/nb-cl-chipseq-results/bam/MundlosLab_IMR-5-75_H3K4me1.trimmed.bwa_hg19.rmdup.filtered.bw", roi=roi, plot.y.min=0, plot.y.max=25, this_color = H3K4me1_color) +
        theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
      make_figure("/Volumes/Elements/nb-cl-chipseq-results/bam/MundlosLab_IMR-5-75_CTCF.trimmed.bwa_hg19.rmdup.filtered.bw", roi=roi, plot.y.min=0, plot.y.max=25, this_color = CTCF_color) +
        theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
      make_figure("/Volumes/Elements/nb-cl-chipseq-results/bam/MundlosLab_IMR-5-75_RAD21.trimmed.bwa_hg19.rmdup.filtered.bw", roi=roi, plot.y.min=0, plot.y.max=15, this_color = RAD21_color) +
        theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
      make_figure("/Volumes/Elements/Virtual4C/IMR5_hg19_canonical_MAPQ30_merged.hic.MYCN4C_Aug28.Normalization_KR_Viewpoint_16075000_16085000_Res_5000.bw", roi=roi, plot.y.min=0, plot.y.max=4, this_color=v4C_color) +
        theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
      make_figure("/Volumes/Elements/MYCNAmplicon/Results/Boeva_lowMYCNExpr_H3K27ac_MeanFC.bw", roi=roi, plot.y.min=0, plot.y.max=15, this_color="gold3") +
        theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
      make_figure("/Volumes/Elements/nb-cl-chipseq-results/bam/Boeva_CLB-GA_GATA3.trimmed.bwa_hg19.rmdup.filtered.bw", roi=roi, plot.y.min=0, plot.y.max=3, this_color=CRC_color) +
        theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
      make_figure("/Volumes/Elements/nb-cl-chipseq-results/bam/Boeva_CLB-GA_PHOX2B.trimmed.bwa_hg19.rmdup.filtered.bw", roi=roi, plot.y.min=0, plot.y.max=5, this_color=CRC_color) +
        theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
      make_figure("/Volumes/Elements/nb-cl-chipseq-results/bam/Boeva_CLB-GA_HAND2.trimmed.bwa_hg19.rmdup.filtered.bw", roi=roi, plot.y.min=0, plot.y.max=5, this_color=CRC_color),
      nrow=13,
      heights = c(0.5,0.5,1,1,1,1,1,1,1,1,1,1,1)),
    height=5, width=2.5, onefile=F, useDingbats=F)
}

# ------------------------------------------------------------------------------
# Plot CHP-212 fragments
# ------------------------------------------------------------------------------

mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl", 
                   mirror = "useast", 
                   GRCh = 37)

window = GRanges(
  seqnames = rep("chr2", 3),
  ranges = IRanges(start = c(11500000, 15000000, 15600000),
                   end   = c(13000000, 15100000, 16300000))
)

for (r in 1:length(window)){
  roi = window[r]
  
  attributes <- c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "gene_biotype")
  filters <- c("chromosome_name","start","end")
  genes = getBM(attributes=attributes, 
                filters=filters, values=list(chromosome=gsub("chr", "", as.character(seqnames(roi))),
                                             start=as.character(start(roi)),
                                             end=as.character(end(roi))), 
                mart=mart) %>% 
    as_tibble() %>% 
    filter(hgnc_symbol != "", gene_biotype == "protein_coding") %>%
    mutate(chromosome_name = paste0("chr", chromosome_name)) %>% 
    dplyr::select(hgnc_symbol, chromosome_name, start_position, end_position)
  colnames(genes) = c("gene", "seqnames", "start", "end")
  genes_df = genes
  genes_df = rbind(genes_df, 
                   data.frame(gene="NULLGENE",seqnames="chr2", start=0, end=1))
  #genes_= makeGRangesFromDataFrame(genes, keep.extra.columns = T)
  genes.fig = genes_df %>%
    ggplot(aes(x=start, y=1)) +
    geom_rect(xmin=genes_df$start, xmax=genes_df$end, ymin=-Inf, ymax=Inf, color=NA, fill="black") + 
    theme_kons2() + 
    theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank()) +
    xlim(start(roi),end(roi))
  
  ggsave(
    paste0("/Volumes/Elements/MYCNAmplicon/Results/CHP212_EpigeneticsFig_Fragments_",
           as.character(seqnames(roi)),"_",as.character(start(roi)),"_", as.character(end(roi)), 
           ".pdf"),
    egg::ggarrange(
      genes.fig +
        theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
      ggplot() + theme_void(),
      make_copynumber_figure("/Volumes/Elements/nb-cl-wgs/nb-cl-wgs-freec/CHP-212/CHP-212_IlluminaWGS.hg19.rmdup.fixRG.bam_CNVs.fillgaps.txt", roi=roi) +
        theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
      make_figure("/Volumes/Elements/nb-cl-atacseq-results/bam/Mundlos_CHP212_ATAC.trimmed.hg19.rmdup.filterednormed.bw", roi=roi, plot.y.min=0, plot.y.max=50, this_color = ATAC_color) +
        theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
      make_figure("/Volumes/Elements/nb-cl-chipseq-results/bam/Boeva_CHP212_H3K27ac.trimmed.bwa_hg19.rmdup.filtered.bw", roi=roi, plot.y.min=0, plot.y.max=30, this_color = H3K27ac_color) +
        theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
      make_figure("/Volumes/Elements/nb-cl-chipseq-results/bam/MundlosLab_CHP-212_H3K4me1.trimmed.bwa_hg19.rmdup.filtered.bw", roi=roi, plot.y.min=0, plot.y.max=15, this_color = H3K4me1_color) +
        theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
      make_figure("/Volumes/Elements/nb-cl-chipseq-results/bam/MundlosLab_CHP-212_CTCF.trimmed.bwa_hg19.rmdup.filtered.bw", roi=roi, plot.y.min=0, plot.y.max=40, this_color = CTCF_color) +
        theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
      make_figure("/Volumes/Elements/Virtual4C/CHP_hg19_canonical_MAPQ30_merged.hic.MYCN4C_Aug28.Normalization_KR_Viewpoint_16075000_16085000_Res_5000.bw", roi=roi, plot.y.min=0, plot.y.max=0.5, this_color=v4C_color) +
        theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
      make_figure("/Volumes/Elements/MYCNAmplicon/Results/Boeva_lowMYCNExpr_H3K27ac_MeanFC.bw", roi=roi, plot.y.min=0, plot.y.max=15, this_color="gold3") +
        theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
      make_figure("/Volumes/Elements/nb-cl-chipseq-results/bam/Boeva_CLB-GA_GATA3.trimmed.bwa_hg19.rmdup.filtered.bw", roi=roi, plot.y.min=0, plot.y.max=3, this_color=CRC_color) +
        theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
      make_figure("/Volumes/Elements/nb-cl-chipseq-results/bam/Boeva_CLB-GA_PHOX2B.trimmed.bwa_hg19.rmdup.filtered.bw", roi=roi, plot.y.min=0, plot.y.max=5, this_color=CRC_color) +
        theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
      make_figure("/Volumes/Elements/nb-cl-chipseq-results/bam/Boeva_CLB-GA_HAND2.trimmed.bwa_hg19.rmdup.filtered.bw", roi=roi, plot.y.min=0, plot.y.max=5, this_color=CRC_color),
      nrow=12,
      heights = c(0.5,0.5,1,1,1,1,1,1,1,1,1,1)),
    height=5, width=2.5, onefile=F, useDingbats=F)
}

# library(Gviz)
# library(GenomicFeatures)
# samplefile <- system.file("extdata", "hg19_knownGene_sample.sqlite", package="GenomicFeatures")
# txdb <- loadDb(samplefile)
# txTr <- GeneRegionTrack(txdb, chromosome="chr2", start=14000000, end=18000000)
# plotTracks(txTr, from=14000000, to=18000000)
# 
# 
# library(ggbio)
# library(Homo.sapiens)
# data(genesymbol, package = "biovizBase")
# wh <- genesymbol[c("MYCN", "NBAS")]
# wh <- range(wh, ignore.strand = TRUE)
# p.txdb <- autoplot(Homo.sapiens, which = wh, stat="reduce")
# p.txdb
# 
# library(EnsDb.Hsapiens.v75)
# ensdb <- EnsDb.Hsapiens.v75
# gr <- GRanges(seqnames = 2, IRanges(14000000, 18000000), strand = "*")
# autoplot(ensdb, GRangesFilter(gr), names.expr = "gene_name", mode="full", stat="reduce") + 
#   theme_kons2()
# 
# genes2.fig

# roi_style = roi
# seqlevelsStyle(roi_style) = "NCBI"
# genes2.fig = 
#   autoplot(ensdb, GRangesFilter(roi_style), names.expr = "gene_name", stat="reduce") + 
#   theme_kons2() + 
#   theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank())
