rm(list=ls())
#devtools::install_local("/Volumes/Elements/MYCNAmplicon/Ext/gTrack-master") # includes my changes
#library(gTrack)
#devtools::load_all("/Volumes/Elements/MYCNAmplicon/Ext/gGnome-master") # includes my changes
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)
library(dplyr)
library(parallel)
library(plyranges)
library(gUtils)
library(gGnome)

# ------------------------------------------------------------------------------
# Read Copy number data
# ------------------------------------------------------------------------------

blacklist_fname = "/Volumes/Elements/MYCNAmplicon/Data/hg19-blacklist.v2.bed"
blacklist = read.table(blacklist_fname, sep="\t", header=F)
colnames(blacklist) = c("Chr", "Start", "End", "Class")
blacklist_gr = makeGRangesFromDataFrame(blacklist, keep.extra.columns = T)

cnv = read.table("/Volumes/Elements/nb-cl-wgs/nb-cl-wgs-freec/IMR-5-75/IMR-5-75_IlluminaWGS.hg19.rmdup.fixRG.bam_CNVsfillgaps.txt",
                 header=F, sep="\t")
colnames(cnv) = c("seqnames", "start", "end", "width", "strand", "copy.number", "copy.number.call")
cnv = makeGRangesFromDataFrame(cnv, keep.extra.columns = T)

wgs_cov = plyranges::read_bigwig("/Volumes/Elements/nb-cl-wgs/nb-cl-wgs-bam/IMR-5-75_IlluminaWGS.hg19.rmdup.fixRG.bw")
wgs_cov_cov = coverage(wgs_cov, weight=wgs_cov$score)

# ------------------------------------------------------------------------------
# Parse svaba manually
# ------------------------------------------------------------------------------

# bcftools index IMR-5-75.svaba.unfiltered.sv.vcf.gz
# bcftools view IMR-5-75.svaba.unfiltered.sv.vcf.gz --regions chr2:10000000-100000000 > IMR-5-75.svaba.unfiltered.sv.vcf.chr2_10M_100M.vcf

# ------------------------------------------------------------------------------
# Create gGraph from junctions 
# ------------------------------------------------------------------------------

out_vcf = "/Volumes/Elements/nb-cl-wgs/playground/IMR-5-75.svaba.unfiltered.sv.vcf.chr2_10M_100M.vcf" # needs to be a vcf

svaba = jJ(rafile=out_vcf, geno = TRUE) # I needed to change utils.R in the gGnome packages (VariantAnnotation::geno instead of geno)

blacklisted_svaba_indices = sort(unique(c(
  queryHits(findOverlaps(svaba$left, blacklist_gr)), 
  queryHits(findOverlaps(svaba$right, blacklist_gr))
  )))
svaba_filtered = svaba[-blacklisted_svaba_indices]

svaba_filtered = svaba_filtered[AD_CELLLINE>50 & SPAN > 1000] 

gg = gG(juncs=svaba_filtered)

# Assign each node with its mean bw signal
gg$nodes$mark(cn=binnedAverage(gg$nodes$gr, wgs_cov_cov, "cn", na.rm=T) %>% .$cn)
gg$nodes$mark(cn_log10=log10(gg$nodes$dt$cn + 0.0001))

gg$set(gr.labelfield = 'node.id')
gg$edges$mark(h=5) # This is how you change the 'curviness'
gg$edges$mark(lwd=0.1*(gg$edges$dt$AD_CELLLINE))
plot(gg$gt, '2:10000000-100000000')

# ------------------------------------------------------------------------------
# Filter out normal copy number segments
# ------------------------------------------------------------------------------
amplificationCutoff = 10 * wgs_cov_cov$chr2 %>% median()
gg = gg[cn>=amplificationCutoff,]

#gg$set(y.field = 'cn')
gg$set(gr.labelfield = 'node.id')
gg$edges$mark(h=1) # This is how you change the 'curviness'
gg$edges$mark(lwd=log10(gg$edges$dt$AD_CELLLINE))
plot(gg$gt, '2:10000000-100000000')

# ------------------------------------------------------------------------------
# Filter out reference edges if ALT edge is of high allele depth
# ------------------------------------------------------------------------------

print(gg$nodes %>% length())
print(gg$edges %>% length())
AD_threshold = floor(quantile(gg$edges$dt$AD_CELLLINE, 0.75, na.rm=T))
#QUAL_threshold = floor(quantile(gg$edges$dt$QUAL, 0.5, na.rm=T))

for (i in 1:length(gg$nodes)){
  print(i)

  this_node = gg$nodes$dt$node.id[i]

  # this_node left
  this_edges = ((gg$edges$dt$n1 == this_node) & (gg$edges$dt$n1.side == "left")) | ((gg$edges$dt$n2 == this_node) & (gg$edges$dt$n2.side == "left"))
  this_edges = gg$edges$dt[this_edges]
  if (any(this_edges$AD_CELLLINE > AD_threshold, na.rm=T)){  # if at this side of the node, there is any deep node
    #gg[,(type == "REF" & (edge.id %in% this_edges$edge.id))]
    gg = gg[,!(type == "REF" & (edge.id %in% this_edges$edge.id))] # filter out all REF edges at this side
  }

  # this_node right
  this_edges = ((gg$edges$dt$n1 == this_node) & (gg$edges$dt$n1.side == "right")) | ((gg$edges$dt$n2 == this_node) & (gg$edges$dt$n2.side == "right"))
  this_edges = gg$edges$dt[this_edges]
  if (any(this_edges$AD_CELLLINE > AD_threshold, na.rm=T)){
    #gg[,(type == "REF" & (edge.id %in% this_edges$edge.id))]
    gg = gg[,!(type == "REF" & (edge.id %in% this_edges$edge.id))]
  }
}

print(gg$nodes %>% length())
print(gg$edges %>% length())

gg$set(gr.labelfield = 'node.id')
gg$edges$mark(h=1) # This is how you change the 'curviness'
gg$edges$mark(lwd=log10(gg$edges$dt$AD_CELLLINE))
plot(gg$gt, '2:10000000-100000000')

# ------------------------------------------------------------------------------
# Create strongly connected clusters
# ------------------------------------------------------------------------------

# # Get strong clusters
# gg$clusters(mode="strong")
# 
# # MYCN containing node
# MYCN_start = 16080683
# MYCN_end = 16087129
# mycn_node = gg$nodes[seqnames == "chr2" & start<= MYCN_start & end>=MYCN_start,]$dt$snode.id
# 
# # MYCN containing cluster
# mycn_cluster = gg$nodes[seqnames == "chr2" & start<= MYCN_start & end>=MYCN_start,]$dt$cluster
# 
# # MYCN containing strong cluster
# gg.sub = gg[cluster==mycn_cluster,]
# 
# # Plot
# gg.sub$set(gr.labelfield = 'node.id')
# gg.sub$edges$mark(h=1) # This is how you change the 'curviness'
# gg.sub$edges$mark(lwd=log10(gg.sub$edges$dt$AD_CELLLINE))
# plot(gg.sub$gt, '2:10000000-100000000')
# 

gg.sub = gg

# Get node ID that contains MYCN (chr2:16,080,683-16,087,129 GRCh37/hg19)
MYCN_start = 16080683
MYCN_end = 16087129
mycn_gg.sub_node = gg.sub$nodes[seqnames == "chr2" & start<= MYCN_start & end>=MYCN_start,]$dt$snode.id

# ------------------------------------------------------------------------------
# Create walks object for neighborhood of MYCN
# ------------------------------------------------------------------------------

#View(gg.sub$nodes$dt)
#View(gg.sub$edges$dt)

# Create walks object
w = gg.sub$walks() # takes 15 min for 22 nodes and 29 edges

# subset of circular walks that contain MYCN
w.mycn = w[sapply(w$snode.id, function (x) (mycn_gg.sub_node %in% x) | (-mycn_gg.sub_node %in% x))] # only walks that contain this_node
w.mycn = w.mycn[circular==TRUE]
w.mycn.nodups = w.mycn[sapply(w.mycn$snode.id, function (sn) sum(duplicated(abs(sn)))==0)]
longest_walk_nodups = w.mycn.nodups[(order(wid, decreasing = T))][1]

# ------------------------------------------------------------------------------
# Overview Plot
# ------------------------------------------------------------------------------

window = GRanges(
  seqnames = rep("chr2", 6),
  ranges = IRanges(start = c(14500000, 17300000, 29500000, 53300000, 66700000, 69100000),
                   end   = c(16200000, 17400000, 30000000, 53700000, 67800000, 70000000))
)
library(BSgenome.Hsapiens.UCSC.hg19)
standardchrs = standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)[1:23]
binsize = 1000 
window_bins = tileGenome(seqlengths = seqlengths(BSgenome.Hsapiens.UCSC.hg19),tilewidth = binsize, cut.last.tile.in.chrom = T) %>%
  filter_by_overlaps(window, minoverlap = binsize/10)
seqlevels(window_bins, pruning.mode="coarse") = standardchrs

Sys.setenv(GENCODE_DIR = "/Volumes/Elements/")
gff = readRDS(gzcon(url('http://mskilab.com/gGnome/hg19/gencode.v19.annotation.gtf.gr.rds')))
genes = gff %Q% (type == 'gene' & gene_name == "MYCN")
seqlevels(genes) = paste0("chr", seqlevels(genes))
gencode = track.gencode(stack.gap = 2e5, cex.label = 0.5, height = 20, name = 'GENCODE', build='hg19', labels.suppress.gr=TRUE, gene.collapse = TRUE, genes = gff %Q% (type == 'gene' & gene_type == "protein_coding") %>% .$gene_name)

gt.copynumber = gTrack(cnv, y.field = "copy.number", y0=0, y1=250)

make_bigwig_gTrack = function(fname, name = "NoName", y0=0, y1=50){
  this_track = read_bigwig(fname, overlap_ranges = window)
  seqlevels(this_track , pruning.mode="coarse") = standardchrs
  this_track = binnedAverage(window_bins, coverage(this_track, weight=this_track$score), "score")
  this_gt = gTrack(this_track, y.field="score", bar=TRUE, name = name, y0=y0, y1=y1) #, 
  return(this_gt)
}

gt.4C = make_bigwig_gTrack("/Volumes/Elements/Virtual4C/IMR5_hg19_canonical_MAPQ30_merged.hic.MYCN4C_Aug28.Normalization_KR_Viewpoint_16075000_16085000_Res_5000.bw", name = "v4C", y0=0, y1=5)
gt.ATAC = make_bigwig_gTrack("/Volumes/Elements/nb-cl-atacseq-results/bam/Mundlos_IMR575_ATAC.trimmed.hg19.rmdup.filterednormed.bw", name = "ATAC", y0=0, y1=50) #, 
gt.H3K27ac = make_bigwig_gTrack("/Volumes/Elements/nb-cl-pe-chipseq-results/bam/Schulte_IMR5_H3K27ac.trimmed.hg19.rmdup.filterednormed.bw", name = "H3K27ac", y0=0, y1=50)
gt.H3K4me1 = make_bigwig_gTrack("/Volumes/Elements/nb-cl-chipseq-results/bam/MundlosLab_IMR-5-75_H3K4me1.trimmed.bwa_hg19.rmdup.filtered.bw", name = "H3K4me1", y0=0, y1=20)
gt.CTCF = make_bigwig_gTrack("/Volumes/Elements/nb-cl-chipseq-results/bam/MundlosLab_IMR-5-75_CTCF.trimmed.bwa_hg19.rmdup.filtered.bw", name = "CTCF", y0=0, y1=20)
plot(c(gt.H3K27ac, gt.CTCF), window)

pdf("/Volumes/Elements/MYCNAmplicon/Results/IMR575_Reconstruction.pdf", height=15, width=10, onefile=F, useDingbats = F)
plot(c(gencode, gt.copynumber, gg.sub$gt, longest_walk_nodups$gt, gt.4C, gt.ATAC, gt.H3K4me1, gt.H3K27ac), window)
dev.off()

# ------------------------------------------------------------------------------
# add hic data
# ------------------------------------------------------------------------------
library(Rcpp)
library(dplyr)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg19)
library(RcppCNPy)
library(data.table)
library(gTrack)
sourceCpp("/Volumes/Transcend/HiC/straw-master/R/straw-R.cpp")

window = GRanges(
  seqnames = rep("chr2", 6),
  ranges = IRanges(start = c(14500000, 17300000, 29500000, 53300000, 66700000, 69100000),
                   end   = c(16200000, 17400000, 30000000, 53700000, 67800000, 70000000))
)

hic_fname = "/Volumes/Elements/MariaSalaDaten/HiC/IMR5/IMR5_hg19_canonical_MAPQ30_merged.hic"
hic_fname = "/Volumes/Elements/MariaSalaDaten/robert_email_20-03-31/IMR_5_75_merged_shallow_sequencing_MAPQ30.hic" # new April

chrs = "chr2" #standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)[1:24] 
binsizekb = 25

buildfullmatrixforchrpair = function(df,binsize,chrR,chrC){
  nr = unname(ceiling(seqlengths(BSgenome.Hsapiens.UCSC.hg19)[chrR] / binsize))
  nc = unname(ceiling(seqlengths(BSgenome.Hsapiens.UCSC.hg19)[chrC] / binsize))
  df[,1] = ((df[,1]) / binsize) + 1
  df[,2] = ((df[,2]) / binsize) + 1
  m = matrix(0L, nrow=nr, ncol=nc, byrow=T)
  if (chrR == chrC){
    for (i in 1:nrow(df)){
      m[df[i,1], df[i,2]] = df[i,3]
      m[df[i,2], df[i,1]] = df[i,3]
    }
  } else {
    for (i in 1:nrow(df)){
      m[df[i,1], df[i,2]] = df[i,3]
    }
  }
  m
}

buildfullmatrixfromhic = function(fname,chrs,binsize,Normalization="NONE"){
  # chrs MUST BE IN ORDER OF STANDARD GENOME (chr1, chr2, ..., chr22, chrX)!!!!! TODO: implement test for that
  # ONLY WORKS for hg19 and allowed binsizes
  bins = ceiling(seqlengths(BSgenome.Hsapiens.UCSC.hg19)[chrs] / binsize)
  bins = unname(bins)
  starti = c(1,cumsum(bins[1:(length(bins)-1)])+1)
  endi = cumsum(bins)
  nr = sum(bins)
  m = matrix(0L,nrow=nr, ncol=nr, byrow=T)
  for (i in 1:length(chrs)){
    for (j in i:length(chrs)){
      m[starti[i]:endi[i],starti[j]:endi[j]] =
        straw_R(paste0(Normalization, " ", fname, " ", gsub("chr", "", chrs[i]), " ", gsub("chr", "", chrs[j]), " BP ", sprintf("%.0f", binsize))) %>%
        buildfullmatrixforchrpair(.,binsize,chrs[i],chrs[j])
      m[starti[j]:endi[j],starti[i]:endi[i]] =
        straw_R(paste0(Normalization, " ", fname, " ", gsub("chr", "", chrs[i]), " ", gsub("chr", "", chrs[j]), " BP ", sprintf("%.0f", binsize))) %>%
        buildfullmatrixforchrpair(.,binsize,chrs[i],chrs[j]) %>%
        t()
    }
  }
  m
}
matrix_labels_granges = function(m, chrs, binsizekb){
  lengths = unname(ceiling(seqlengths(BSgenome.Hsapiens.UCSC.hg19)[chrs] / (binsizekb*1000)))
  labels_chr = NA*1:ncol(m)
  labels_start = NA*1:ncol(m)
  labels_end = NA*1:ncol(m)
  for (ci1 in 1:length(chrs)){
    start_c1 = sum(lengths[0:(ci1-1)])+1
    end_c1 = sum(lengths[0:ci1])
    labels_chr[start_c1:end_c1] = chrs[ci1]
    labels_start[start_c1:end_c1] = (0:(lengths[ci1]-1))*binsizekb*1000
    labels_end[start_c1:end_c1] = ((0:(lengths[ci1]-1))+1)*binsizekb*1000-1
  }
  lab = GRanges(seqnames=labels_chr, ranges=IRanges(start = labels_start, end=labels_end), strand = "*")
  return(lab)
}
m_kr = buildfullmatrixfromhic(fname = hic_fname,
                              chrs = chrs,
                              binsize = binsizekb*1000, 
                              Normalization = "KR")
m_kr_labels = matrix_labels_granges(m_kr, chrs, binsizekb)
gt.hic = gTrack(m_kr_labels, mdata = m_kr, colormap = c('white', 'red', 'black'), height=100)


# Create gTrack for all non-blacklisted rearrangements
svaba_forplotting = jJ(rafile=out_vcf, geno = TRUE) # I needed to change utils.R in the gGnome packages (VariantAnnotation::geno instead of geno)
blacklisted_svaba_forplotting_indices = sort(unique(c(
  queryHits(findOverlaps(svaba_forplotting$left, blacklist_gr)), 
  queryHits(findOverlaps(svaba_forplotting$right, blacklist_gr))
)))
svaba_forplotting = svaba[-blacklisted_svaba_forplotting_indices]
gg_forplotting = gG(juncs=svaba_forplotting)
gg_forplotting$edges$mark(h=3) # This is how you change the 'curviness'
gg_forplotting$edges$mark(lwd=0.05*(gg_forplotting$edges$dt$AD_CELLLINE))
plot(sv.gt, window)
sv.gt = gg_forplotting$gt
sv.gt$height = 30

# Create gTrack for reconstruction
reconstruction.gt = longest_walk_nodups$gt
reconstruction.gt$height = 20

# Load the additional tracks
gt.NoMYCNExpr_H3K27ac = make_bigwig_gTrack("/Volumes/Elements/MYCNAmplicon/Results/Boeva_nonMYCNExpr_H3K27ac_MeanFC.bw", name = "Non MYCN Expr NB (H3K27ac FC)", y0=0, y1=20)
gt.LowMYCNExpr_H3K27ac = make_bigwig_gTrack("/Volumes/Elements/MYCNAmplicon/Results/Boeva_lowMYCNExpr_H3K27ac_MeanFC.bw", name = "Low MYCN Expr NB (H3K27ac FC)", y0=0, y1=20)
gt.CLBGA_GATA3 = make_bigwig_gTrack("/Volumes/Elements/nb-cl-chipseq-results/bam/Boeva_CLB-GA_GATA3.trimmed.bwa_hg19.rmdup.filtered.bw", name = "CLB-GA GATA3", y0=0, y1=10)
gt.CLBGA_PHOX2B = make_bigwig_gTrack("/Volumes/Elements/nb-cl-chipseq-results/bam/Boeva_CLB-GA_PHOX2B.trimmed.bwa_hg19.rmdup.filtered.bw", name = "CLB-GA PHOX2B", y0=0, y1=10)
gt.CLBGA_HAND2 = make_bigwig_gTrack("/Volumes/Elements/nb-cl-chipseq-results/bam/Boeva_CLB-GA_HAND2.trimmed.bwa_hg19.rmdup.filtered.bw", name = "CLB-GA HAND2", y0=0, y1=10)

pdf("/Volumes/Elements/MYCNAmplicon/Results/IMR575_Reconstruction2_April6.pdf", height=14, width=7, onefile=F, useDingbats = F)
plot(c(gencode, gt.LowMYCNExpr_H3K27ac, gt.copynumber, sv.gt, reconstruction.gt, gt.ATAC, gt.H3K27ac, gt.CTCF, gt.4C, gt.hic), window,
     y.grid.col=F, y.grid.lwd=F, sep.lwd=0, sep.col=F, yaxis.pretty = 1)
dev.off()
#save.image("/Volumes/Elements/MYCNAmplicon/Workspace/IMR575Reconstruction.Rdata")
#load("/Volumes/Elements/MYCNAmplicon/Workspace/IMR575Reconstruction.Rdata")

save.image("/Volumes/Elements/MYCNAmplicon/Workspace/IMR575Reconstruction_April6.Rdata")
load("/Volumes/Elements/MYCNAmplicon/Workspace/IMR575Reconstruction_April6.Rdata")


# ------------------------------------------------------------------------------
se = read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_lowMYCN_CRCdriven_SE.bed")
strand(se) = "*"
gt.se = gTrack(se, col="gold2", height=5)
reconstruction.gt = longest_walk_nodups$gt
reconstruction.gt$height = 15
gencode = track.gencode(stack.gap = 2e5, cex.label = 0.5, height = 5, name = 'GENCODE', build='hg19', labels.suppress.gr=TRUE, gene.collapse = TRUE, genes = gff %Q% (type == 'gene' & gene_type == "protein_coding") %>% .$gene_name)
make_bigwig_gTrack = function(fname, name = "NoName", y0=0, y1=50){
  this_track = read_bigwig(fname, overlap_ranges = window)
  seqlevels(this_track , pruning.mode="coarse") = standardchrs
  this_track = binnedAverage(window_bins, coverage(this_track, weight=this_track$score), "score")
  this_gt = gTrack(this_track, y.field="score", bar=TRUE, name = name, y0=y0, y1=y1) #, 
  return(this_gt)
}
gt.CTCF = make_bigwig_gTrack("/Volumes/Elements/nb-cl-chipseq-results/bam/MundlosLab_IMR-5-75_CTCF.trimmed.bwa_hg19.rmdup.filtered.bw", name = "CTCF", y0=0, y1=20)
#pdf("/Volumes/Elements/MYCNAmplicon/Results/IMR575_Reconstruction_Dec11.pdf", height=12, width=7, onefile=F, useDingbats = F)
pdf("/Volumes/Elements/MYCNAmplicon/Results/IMR575_Reconstruction_FINAL_April6.pdf", height=12, width=7, onefile=F, useDingbats = F)
plot(c(gencode, gt.se, gt.copynumber, reconstruction.gt, gt.H3K27ac, gt.CTCF, gt.4C, gt.hic), window,
     y.grid.col=F, y.grid.lwd=F, sep.lwd=0, sep.col=F, yaxis.pretty = 1)
dev.off()

# longest_walk_nodups$grl %>% 
#   as_tibble() %>% 
#   write.table("/Volumes/Elements/MYCNAmplicon/Results/IMR575_Reconstruction.txt",
#               sep="\t", quote=F)
