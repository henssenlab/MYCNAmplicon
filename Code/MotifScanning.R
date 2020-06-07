library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(JASPAR2018) 
library(TFBSTools)
library(plyranges)
library(dplyr)
library(ggplot2)
library(egg)
library(RColorBrewer)
source("/Volumes/Elements/MYCNAmplicon/Code/CustomThemes.R")

getMotifHits = function(in_bed_file, 
                        out_bed_file, 
                        tf_name, 
                        filter_chr=standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)[1:23], 
                        filter_start=-Inf,
                        filter_end=Inf, 
                        margin = 0,
                        species_id = 9606){
  
  # in_bed_file = "/Volumes/Elements/MYCNAmplicon/Results/Boeva_nMNA_MYCN_Enhancers.bed"
  # out_bed_file = NULL
  # tf_name = "MYCN"
  # filter_chr=standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)[1:23]
  # filter_start=-Inf
  # filter_end=Inf
  # margin = 0
  # tf_name = "Isl1"
  # filter_chr=standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)[1:23]
  # species_id = 10090
  
  intervals = read_bed(as.character(in_bed_file)) %>%
    filter(seqnames %in% filter_chr,
           start >= filter_start, 
           end <= filter_end)
  start(intervals) = start(intervals) - margin
  end(intervals) = end(intervals) + margin
  
  seqlevels(intervals, pruning.mode="coarse") = seqlevels(BSgenome.Hsapiens.UCSC.hg19)
  seqinfo(intervals) = seqinfo(BSgenome.Hsapiens.UCSC.hg19)
  seq = getSeq(BSgenome.Hsapiens.UCSC.hg19, intervals, as.character=F)
  
  pfm <- getMatrixSet(JASPAR2018, list(species=species_id, name=tf_name))[[1]]
  
  hits_fwd <- lapply(seq, function(s) matchPWM(as.matrix(pfm), s, min.score="85%", with.score = TRUE) )
  hits_rev <- lapply(seq, function(s) matchPWM(reverseComplement(as.matrix(pfm)), s, min.score="85%", with.score = TRUE) )
  
  motif_hits = GRanges()
  for (i in 1:length(hits_fwd)){
    
    isForward = TRUE
    
    this_starts = start(hits_fwd[[i]])
    this_ends = end(hits_fwd[[i]])
    this_scores = mcols(hits_fwd[[i]])$score
    
    if (length(this_starts) == 0) next
    if (length(this_starts) != length(this_ends)) stop("vector this_starts must be the same length as this_ends")
    
    this_chr = rep(as.character(seqnames(intervals)[i]), length(this_starts))
    this_starts = this_starts + start(intervals)[i] -2  # - 1 would make more sense, but -2 matches IGV
    this_ends = this_ends + start(intervals)[i] - 2 # - 1
    this_strand = rep(ifelse(isForward, "+", "-"), length(this_starts))
    
    motif_hits = c(motif_hits, 
                   GRanges(seqnames = this_chr,
                           ranges = IRanges(start = this_starts, end = this_ends),
                           strand = this_strand, 
                           mcols = data.frame("score" = this_scores)))
    
  }
  
  for (i in 1:length(hits_rev)){
    
    isForward = FALSE
    
    this_starts = start(hits_rev[[i]])
    this_ends = end(hits_rev[[i]])
    this_scores = mcols(hits_rev[[i]])$score
    
    if (length(this_starts)==0) next
    if (length(this_starts) != length(this_ends)) stop("vector this_starts must be the same length as this_ends")
    
    this_chr = rep(as.character(seqnames(intervals)[i]), length(this_starts))
    this_starts = this_starts + start(intervals)[i] - 1
    this_ends = this_ends + start(intervals)[i] - 1
    this_strand = rep(ifelse(isForward, "+", "-"), length(this_starts))
    
    motif_hits = c(motif_hits, 
                   GRanges(seqnames = this_chr,
                           ranges = IRanges(start = this_starts, end = this_ends),
                           strand = this_strand,
                           mcols = data.frame("score" = this_scores)))
  }
  
  if (length(motif_hits)>0) motif_hits$TF = tf_name
  
  #write_bed(x=motif_hits %>% filter(strand == "+"), file=gsub(".bed$", "fwd.bed", out_bed_file))
  #write_bed(x=motif_hits %>% filter(strand == "-"), file=gsub(".bed$", "rev.bed", out_bed_file))
  #write_bed(x = motif_hits, file = out_bed_file)
  
  return(motif_hits)
}

getMotifHitsFromMotifFile = function(in_bed_file, 
                                     out_bed_file, 
                                     tf_name,
                                     pfm_file, 
                                     filter_chr=standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)[1:23], 
                                     filter_start=-Inf,
                                     filter_end=Inf, 
                                     margin = 0){
  
  # in_bed_file = "/Volumes/Elements/MYCNAmplicon/Results/Boeva_nMNA_MYCN_Enhancers.bed"
  # out_bed_file = NULL
  # tf_name = "MYCN"
  # filter_chr=standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)[1:23]
  # filter_start=-Inf
  # filter_end=Inf
  # margin = 0
  # tf_name = "Isl1"
  # pfm_file ="/Volumes/Elements/MYCNAmplicon/Data/MA1608.1.jaspar" 
  
  intervals = read_bed(as.character(in_bed_file)) %>%
    filter(seqnames %in% filter_chr,
           start >= filter_start, 
           end <= filter_end)
  start(intervals) = start(intervals) - margin
  end(intervals) = end(intervals) + margin
  seqlevels(intervals, pruning.mode="coarse") = seqlevels(BSgenome.Hsapiens.UCSC.hg19)
  seqinfo(intervals) = seqinfo(BSgenome.Hsapiens.UCSC.hg19)
  seq = getSeq(BSgenome.Hsapiens.UCSC.hg19, intervals, as.character=F)
  pfm <- readJASPARMatrix(pfm_file)
  hits_fwd <- lapply(seq, function(s) matchPWM(as.matrix(pfm), s, min.score="85%", with.score = TRUE) )
  hits_rev <- lapply(seq, function(s) matchPWM(reverseComplement(as.matrix(pfm)), s, min.score="85%", with.score = TRUE) )
  motif_hits = GRanges()
  for (i in 1:length(hits_fwd)){
    isForward = TRUE
    this_starts = start(hits_fwd[[i]])
    this_ends = end(hits_fwd[[i]])
    this_scores = mcols(hits_fwd[[i]])$score
    if (length(this_starts) == 0) next
    if (length(this_starts) != length(this_ends)) stop("vector this_starts must be the same length as this_ends")
    this_chr = rep(as.character(seqnames(intervals)[i]), length(this_starts))
    this_starts = this_starts + start(intervals)[i] -2  # - 1 would make more sense, but -2 matches IGV
    this_ends = this_ends + start(intervals)[i] - 2 # - 1
    this_strand = rep(ifelse(isForward, "+", "-"), length(this_starts))
    motif_hits = c(motif_hits, 
                   GRanges(seqnames = this_chr,
                           ranges = IRanges(start = this_starts, end = this_ends),
                           strand = this_strand, 
                           mcols = data.frame("score" = this_scores)))
  }
  for (i in 1:length(hits_rev)){
    isForward = FALSE
    this_starts = start(hits_rev[[i]])
    this_ends = end(hits_rev[[i]])
    this_scores = mcols(hits_rev[[i]])$score
    if (length(this_starts)==0) next
    if (length(this_starts) != length(this_ends)) stop("vector this_starts must be the same length as this_ends")
    this_chr = rep(as.character(seqnames(intervals)[i]), length(this_starts))
    this_starts = this_starts + start(intervals)[i] - 1
    this_ends = this_ends + start(intervals)[i] - 1
    this_strand = rep(ifelse(isForward, "+", "-"), length(this_starts))
    motif_hits = c(motif_hits, 
                   GRanges(seqnames = this_chr,
                           ranges = IRanges(start = this_starts, end = this_ends),
                           strand = this_strand,
                           mcols = data.frame("score" = this_scores)))
  }
  if (length(motif_hits)>0) motif_hits$TF = tf_name
  return(motif_hits)
}

# ------------------------------------------------------------------------------
# Search regions of interest for all TF in JASPAR2018
# ------------------------------------------------------------------------------

tf = lapply(getMatrixSet(JASPAR2018, list(species=9606)), name) %>% unname() %>% unlist() # a list of all names of TFs in the JASPAR2018 database
hits = list()
for (this_tf in tf){
  print(this_tf)
  hits[[this_tf]] = getMotifHits("/Volumes/Elements/MYCNAmplicon/Results/Boeva_nMNA_MYCN_Enhancers.bed",
                                 NULL,
                                 this_tf)
}
hits_gr = do.call(c, unname(hits))
hits_gr_all = hits_gr

enhancers = read_bed(as.character("/Volumes/Elements/MYCNAmplicon/Results/Boeva_nMNA_MYCN_Enhancers.bed")) 
hits_gr = hits_gr %>%
  join_overlap_inner(enhancers %>% plyranges::select(name)) 
hits_gr %>% 
  as_tibble() %>% 
  group_by(TF) %>%
  summarise(n_hits = dplyr::n(),
            n_regions = n_distinct(name),
            regions = paste(unique(name), collapse = ",")) %>%
  write.table("/Volumes/Elements/MYCNAmplicon/Results/EnhancerRegionsMotifHits.txt",
              quote=F, sep="\t", row.names = F, col.names = T)


# ------------------------------------------------------------------------------
# ISL1
# ------------------------------------------------------------------------------

isl1_hits = getMotifHitsFromMotifFile("/Volumes/Elements/MYCNAmplicon/Results/Boeva_nMNA_MYCN_Enhancers.bed",
                                      NULL,
                                      tf_name = "Isl1",
                                      pfm_file = "/Volumes/Elements/MYCNAmplicon/Data/JASPAR2020_Isl1_MA1608.1.jaspar",
                                      filter_chr=standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)[1:23], 
                                      filter_start=-Inf,
                                      filter_end=Inf, 
                                      margin = 0)
enhancers = read_bed(as.character("/Volumes/Elements/MYCNAmplicon/Results/Boeva_nMNA_MYCN_Enhancers.bed")) 
isl1_hits = isl1_hits %>%
  join_overlap_inner(enhancers %>% plyranges::select(name)) 
isl1_hits %>% 
  as_tibble() %>% 
  group_by(TF) %>%
  summarise(n_hits = dplyr::n(),
            n_regions = n_distinct(name),
            regions = paste(unique(name), collapse = ",")) %>%
  write.table("/Volumes/Elements/MYCNAmplicon/Results/EnhancerRegionsIsl1Hits.txt",
              quote=F, sep="\t", row.names = F, col.names = T)

# ------------------------------------------------------------------------------
# PHOX2B
# ------------------------------------------------------------------------------

phox2b_hits = getMotifHitsFromMotifFile("/Volumes/Elements/MYCNAmplicon/Results/Boeva_nMNA_MYCN_Enhancers.bed",
                                        NULL,
                                        tf_name = "PHOX2B",
                                        pfm_file = "/Volumes/Elements/MYCNAmplicon/Data/JASPAR2020_PHOX2B_MA0681.2.jaspar",
                                        filter_chr=standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)[1:23], 
                                        filter_start=-Inf,
                                        filter_end=Inf, 
                                        margin = 0)
enhancers = read_bed(as.character("/Volumes/Elements/MYCNAmplicon/Results/Boeva_nMNA_MYCN_Enhancers.bed")) 
phox2b_hits = phox2b_hits %>%
  join_overlap_inner(enhancers %>% plyranges::select(name)) 
phox2b_hits %>% 
  as_tibble() %>% 
  group_by(TF) %>%
  summarise(n_hits = dplyr::n(),
            n_regions = n_distinct(name),
            regions = paste(unique(name), collapse = ",")) %>%
  write.table("/Volumes/Elements/MYCNAmplicon/Results/EnhancerRegionsPHOX2BHits.txt",
              quote=F, sep="\t", row.names = F, col.names = T)

# ------------------------------------------------------------------------------
# HAND2
# ------------------------------------------------------------------------------

hand2_hits = getMotifHitsFromMotifFile("/Volumes/Elements/MYCNAmplicon/Results/Boeva_nMNA_MYCN_Enhancers.bed",
                                       NULL,
                                       tf_name = "HAND2",
                                       pfm_file = "/Volumes/Elements/MYCNAmplicon/Data/JASPAR2020_HAND2_MA1638.1.jaspar",
                                       filter_chr=standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)[1:23], 
                                       filter_start=-Inf,
                                       filter_end=Inf, 
                                       margin = 0)
enhancers = read_bed(as.character("/Volumes/Elements/MYCNAmplicon/Results/Boeva_nMNA_MYCN_Enhancers.bed")) 
hand2_hits = hand2_hits %>%
  join_overlap_inner(enhancers %>% plyranges::select(name)) 
hand2_hits %>% 
  as_tibble() %>% 
  group_by(TF) %>%
  summarise(n_hits = dplyr::n(),
            n_regions = n_distinct(name),
            regions = paste(unique(name), collapse = ",")) %>%
  write.table("/Volumes/Elements/MYCNAmplicon/Results/EnhancerRegionsHAND2Hits.txt",
              quote=F, sep="\t", row.names = F, col.names = T)



# ------------------------------------------------------------------------------
# Search for E-boxes
# ------------------------------------------------------------------------------

intervals = read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_nMNA_MYCN_Enhancers.bed") 
seqlevels(intervals, pruning.mode="coarse") = seqlevels(BSgenome.Hsapiens.UCSC.hg19)
seqinfo(intervals) = seqinfo(BSgenome.Hsapiens.UCSC.hg19)
seq = getSeq(BSgenome.Hsapiens.UCSC.hg19, intervals, as.character=F)
names(seq) = intervals$name
canonicalEBox = vmatchPattern("CACGTG", seq)
noncanonicalEBox = vmatchPattern("CANNTG", seq, fixed=FALSE)

# How many hits in each region
lapply(canonicalEBox, length)
lapply(noncanonicalEBox, length)

# How many hits in each region per kb
sapply(canonicalEBox, length) / (width(intervals) / 1000)
sapply(noncanonicalEBox, length) / (width(intervals) / 1000)

# What is the background for chr2?
mycn_chr = "chr2"
mycn_start = 16080683
mycn_end = 16087129
window_of_interest = 500000
mycn_neighborhood = GRanges(
  seqnames = mycn_chr,
  ranges = IRanges(start = 1,
                   end = seqlengths(BSgenome.Hsapiens.UCSC.hg19)[2]),
  strand = "*"
)
seqlevels(mycn_neighborhood, pruning.mode="coarse") = seqlevels(BSgenome.Hsapiens.UCSC.hg19)
seqinfo(mycn_neighborhood) = seqinfo(BSgenome.Hsapiens.UCSC.hg19)
mycn_neighborhood_seq = getSeq(BSgenome.Hsapiens.UCSC.hg19, mycn_neighborhood, as.character=F)
mycn_neighborhood_canonicalEBox = vmatchPattern("CACGTG", mycn_neighborhood_seq)
mycn_neighborhood_noncanonicalEBox = vmatchPattern("CANNTG", mycn_neighborhood_seq, fixed=FALSE)
sapply(mycn_neighborhood_canonicalEBox, length) / (width(mycn_neighborhood) / 1000)
sapply(mycn_neighborhood_noncanonicalEBox, length) / (width(mycn_neighborhood) / 1000)

# ------------------------------------------------------------------------------
# Plot TF binding sites
# ------------------------------------------------------------------------------
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
  xlim(start(roi), end(roi))
cRE = read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_nMNA_MYCN_Enhancers.bed") %>% as_tibble()
cRE.fig = cRE %>%
  ggplot(aes(x=start, y=score)) +
  geom_rect(aes(xmin=start, xmax=end), ymin=-Inf, ymax=Inf, color="firebrick3", fill="firebrick3") + 
  theme_kons2() + 
  theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank()) +
  xlim(start(roi), end(roi))
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
  xlim(start(roi), end(roi))




plot_tf = function(gr, tf_name, roi, tf_color="black") {
  fig = 
    gr %>%
    filter(seqnames == seqnames(roi)) %>% 
    filter(TF == tf_name) %>% 
    as_tibble() %>%
    ggplot(aes(xmin=start, xmax = start, ymin=0, ymax=1)) + 
    geom_rect(color=tf_color,fill=tf_color, size=0.1) +
    guides(color=F) +
    theme_kons2() +
    xlim(start(roi), end(roi)) +
    xlab(seqnames(roi)) +
    theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank())
  return(fig)
}

ggsave("/Volumes/Elements/MYCNAmplicon/Results/TFMotifs.pdf",
       egg::ggarrange(
         genes.fig + theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
         ggplot() + theme_void(),
         cRE.fig + theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
         ggplot() + theme_void(),
         plot_tf(phox2b_hits, "PHOX2B", roi, tf_color=RColorBrewer::brewer.pal(8, "Spectral")[1]) + theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
         plot_tf(hits_gr_all, "GATA3", roi, tf_color=RColorBrewer::brewer.pal(8, "Spectral")[2]) + theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
         plot_tf(hand2_hits, "HAND2", roi, tf_color=RColorBrewer::brewer.pal(8, "Spectral")[3]) + theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()),
         plot_tf(isl1_hits, "Isl1", roi, tf_color=RColorBrewer::brewer.pal(8, "Spectral")[4]) + theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
         plot_tf(hits_gr_all, "TBX2", roi, tf_color=RColorBrewer::brewer.pal(8, "Spectral")[5]) + theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
         plot_tf(hits_gr_all, "ASCL1", roi, tf_color=RColorBrewer::brewer.pal(8, "Spectral")[6]) + theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
         plot_tf(hits_gr_all, "MYCN", roi, tf_color=RColorBrewer::brewer.pal(8, "Spectral")[7]) + theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
         plot_tf(hits_gr_all, "TEAD4", roi, tf_color=RColorBrewer::brewer.pal(8, "Spectral")[8]),
         nrow = 12,
         heights = c(1,0.5,1,0.25, rep(1,8)),
         labels = c("Genes","", "Enh", "", "PHOX2B", "GATA3", "HAND2", "ISL1", "TBX2", "ASCL1", "MYCN", "TEAD4"),
         label.args = list(gp=grid::gpar(fontface="italic", fontfamily="Helvetica", fontsize=6))
       ),
       height = 3, width=3, 
       useDingbats=F, onefile=F)

# ------------------------------------------------------------------------------
# Overlap with TF ChIP-seq
# ------------------------------------------------------------------------------
library(dplyr)
library(plyranges)

phox2b_peaks = read_narrowpeaks("/Volumes/Elements/nb-cl-chipseq-results/MACS2/Boeva_CLB-GA_PHOX2B/Boeva_CLB-GA_PHOX2B_MACS2_nocontrol_peaks.narrowPeak")
gata3_peaks = read_narrowpeaks("/Volumes/Elements/nb-cl-chipseq-results/MACS2/Boeva_CLB-GA_GATA3/Boeva_CLB-GA_GATA3_MACS2_nocontrol_peaks.narrowPeak")
hand2_peaks = read_narrowpeaks("/Volumes/Elements/nb-cl-chipseq-results/MACS2/Boeva_CLB-GA_HAND2/Boeva_CLB-GA_HAND2_MACS2_nocontrol_peaks.narrowPeak")

enhancers = read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_nMNA_MYCN_Enhancers.bed")

enhancers$Overlaps_PHOX2B_Peak = FALSE
enhancers[unique(queryHits(findOverlaps(enhancers, phox2b_peaks)))]$Overlaps_PHOX2B_Peak = TRUE

enhancers$Overlaps_HAND2_Peak = FALSE
enhancers[unique(queryHits(findOverlaps(enhancers, hand2_peaks)))]$Overlaps_HAND2_Peak = TRUE

enhancers$Overlaps_GATA3_Peak = FALSE
enhancers[unique(queryHits(findOverlaps(enhancers, gata3_peaks)))]$Overlaps_GATA3_Peak = TRUE

print(enhancers)

enhancers %>% as_tibble() %>% write.table("/Volumes/Elements/MYCNAmplicon/Results/EnhancersVsCRCChIPseqPeaks.txt",
                                          quote=F, sep="\t", row.names=F, col.names=T)


