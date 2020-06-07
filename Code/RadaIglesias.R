library(dplyr)
library(plyranges)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg19)
library(biomaRt)
source("/Volumes/Elements/MYCNAmplicon/Code/CustomThemes.R")
source("/Volumes/Elements/MYCNAmplicon/Code/hasOverlap.R")

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

roi = GRanges(
  seqnames = "chr2",
  ranges = IRanges(start = 16080683-500000,
                   end = 16087129+500000),
  strand = "*"
)

# Genes Figure
mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl", 
                   mirror = "useast", 
                   GRCh = 37)
attributes <- c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "gene_biotype")
filters <- c("chromosome_name","start","end")
genes = getBM(attributes=attributes, 
              filters=filters, values=list(chromosome=gsub("chr", "", as.character(seqnames(roi))),
                                           start=as.character(start(roi)),
                                           end=as.character(end(roi))), 
              mart=mart) %>% 
  as_tibble() %>% 
  filter(hgnc_symbol != "") %>%
  #filter(hgnc_symbol != "", gene_biotype == "protein_coding") %>%
  mutate(chromosome_name = paste0("chr", chromosome_name)) %>% 
  dplyr::select(hgnc_symbol, chromosome_name, start_position, end_position)
colnames(genes) = c("gene", "seqnames", "start", "end")
genes_df = genes
genes_df = rbind(genes_df, 
                 data.frame(gene="NULLGENE",seqnames="chr2", start=0, end=1))
genes_df = 
  data.frame(
    "start" = c(15307032, 15731302, 16060521, 16080683, 16190549, 16730727),
    "end" = c(15701454, 15771235, 16076139, 16087129, 16225923, 16847599),
    "name" = c("NBAS", "DDX1", "MYCNUT", "MYCN", "GACAT3", "FAM49A"),
    "class" = rep("gene", 6)
  )
genes.fig = genes_df %>%
  ggplot(aes(x=start, y=1)) +
  geom_rect(xmin=genes_df$start, xmax=genes_df$end, ymin=-Inf, ymax=Inf, color="black", fill="black", size=0.2) + 
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

ggsave(
  "/Volumes/Elements/MYCNAmplicon/Results/RadaIglesiasH3K27ac_Feb28.pdf",
  egg::ggarrange(
    genes.fig +
      theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
    ggplot() + theme_void(),
    cRE.fig +
      theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
    ggplot() + theme_void(),
    make_figure("/Volumes/Elements/MYCNAmplicon/Results/Boeva_lowMinusNon_H3K27ac_DeltaMeanFC.bw", roi=roi, plot.y.min=0, plot.y.max=15, this_color = "black") +
      theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
    make_figure("/Volumes/Elements/nb-cl-chipseq-results/bam/Rada_ESC_H3K27ac.trimmed.bwa_hg19.rmdup.filtered.bw", roi=roi, plot.y.min=0, plot.y.max=5, this_color = "steelblue") +
      theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
    make_figure("/Volumes/Elements/nb-cl-chipseq-results/bam/Rada_NEC_H3K27ac.trimmed.bwa_hg19.rmdup.filtered.bw", roi=roi, plot.y.min=0, plot.y.max=5, this_color = "steelblue") +
      theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
    make_figure("/Volumes/Elements/nb-cl-chipseq-results/bam/Rada_NCC_H3K27ac.trimmed.bwa_hg19.rmdup.filtered.bw", roi=roi, plot.y.min=0, plot.y.max=5, this_color="steelblue") +
      theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
    make_figure("/Volumes/Elements/nb-cl-chipseq-results/bam/Roadmap_FetalAdrenal_H3K27ac.trimmed.bwa_hg19.rmdup.filtered.bw", roi=roi, plot.y.min=0, plot.y.max=5, this_color="steelblue"),
    #  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
    #make_figure("/Volumes/Elements/nb-cl-chipseq-results/bam/Boeva_hNCC_rep1_H3K27ac.trimmed.bwa_hg19.rmdup.filtered.bw", roi=roi, plot.y.min=0, plot.y.max=5, this_color="steelblue") +
    #  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
    #make_figure("/Volumes/Elements/nb-cl-chipseq-results/bam/Boeva_hNCC_rep2_H3K27ac.trimmed.bwa_hg19.rmdup.filtered.bw", roi=roi, plot.y.min=0, plot.y.max=5, this_color="steelblue"),
    nrow=9,
    heights = c(0.25,0.25,0.25,0.25,1,1,1,1,1)),
  height=2.2, width=2.5, onefile=F, useDingbats=F)
 
