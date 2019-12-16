library(gTrack)
devtools::load_all("/Volumes/Elements/nb-reconstruct-amplicon/scripts/gGnome-master") # includes my changes
#devtools::load_all("/Volumes/Elements/nb-reconstruct-amplicon/scripts/gTrack-master") # includes my changes
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)
library(dplyr)
library(parallel)
library(plyranges)
cores = detectCores()

# ------------------------------------------------------------------------------
# Read Copy number data
# ------------------------------------------------------------------------------

blacklist_fname = "~/Downloads/hg19-blacklist.v2.bed"
#blacklist_fname = "/fast/users/helmsauk_c/work/resources/blacklists/hg19-blacklist.v2.bed"
blacklist = read.table(blacklist_fname, sep="\t", header=F)
colnames(blacklist) = c("Chr", "Start", "End", "Class")
blacklist$Chr = gsub("chr", "", blacklist$Chr)
blacklist_gr = makeGRangesFromDataFrame(blacklist, keep.extra.columns = T)

load("/Volumes/Transcend/PalmTrees/Analysis/WorkspaceData/ascat_cnv.Rdata")
seqlevelsStyle(ascat_cnv_gr) <- "NCBI" # get rid of chr in seqlevels
wgs_cov_cov = ascat_cnv_gr %>% filter(Sample == "NBL26") %>%
  coverage(weight = "TumorTotalCopyNumber")

# ------------------------------------------------------------------------------
# Prepare svaba vcf
# ------------------------------------------------------------------------------

# Filter vcf files 
# changed "ID=PL,Number=." to "ID=PL,Number=G" in header
# bgzip somatic_run.svaba.unfiltered.somatic.sv.vcf
# bcftools index somatic_run.svaba.unfiltered.somatic.sv.vcf.gz
# bcftools view somatic_run.svaba.unfiltered.somatic.sv.vcf.gz --regions 2:10000000-75000000 > NBL26.svaba.unfiltered.somatic.sv.vcf.chr2_10M_75M.vcf

# ------------------------------------------------------------------------------
# Create gGraph from junctions 
# ------------------------------------------------------------------------------

out_vcf = "/Volumes/Elements/MYCNAmplicon/Data/SvABA/NBL26/NBL26.svaba.unfiltered.somatic.sv.vcf.chr2_10M_75M.vcf" # needs to be a vcf

svaba = jJ(rafile=out_vcf, geno = TRUE) # I need to change utils.R in the gGnome packages (VariantAnnotation::geno instead of geno)
blacklisted_svaba_indices = sort(unique(c(
  queryHits(findOverlaps(svaba$left, blacklist_gr)), # those are the filtered bins
  queryHits(findOverlaps(svaba$right, blacklist_gr))
)))
svaba_filtered = svaba[-blacklisted_svaba_indices]

median_AD = median(svaba_filtered$dt$AD_TUMOR, na.rm=T)
#svaba_filtered = svaba[FILTER == "PASS"]
svaba_filtered = svaba_filtered[AD_TUMOR>10*median_AD & SPAN > 1000] 

gg = gG(juncs=svaba_filtered)

# Assign each node with its mean copy number signal
load("/Volumes/Transcend/PalmTrees/Analysis/WorkspaceData/ascat_cnv.Rdata")
seqlevelsStyle(ascat_cnv_gr) <- "NCBI" # get rid of chr in seqlevels
seqlevels(ascat_cnv_gr, pruning.mode="coarse") = seqlevels(gg$nodes$gr)
seqlengths(ascat_cnv_gr) = seqlengths(gg$nodes$gr)
wgs_cov_cov = ascat_cnv_gr %>% filter(Sample == "NBL26") %>%
  coverage(weight = "TumorTotalCopyNumber")
gg$nodes$mark(cn=binnedAverage(gg$nodes$gr, wgs_cov_cov, "cn", na.rm=T) %>% .$cn)
gg$nodes$mark(cn_log10=log10(gg$nodes$dt$cn + 0.0001))

# #lapply(seq_along(gg$nodes$gr), function(i) wgs_cov %>% filter_by_overlaps(gg$nodes$gr[i]) %>% .$score %>% median)
# 
# #gg$nodes$dt %>% View
# #gg$set(y.field = 'cn_log10')
# gg$set(gr.labelfield = 'node.id')
# gg$edges$mark(h=5) # This is how you change the 'curviness'
# gg$edges$mark(lwd=0.1*(gg$edges$dt$AD_TUMOR))
# #plot(gg$gt, '2:15000000-17000000')
# plot(gg$gt, '2:10000000-100000000')
# 
# 
# # ------------------------------------------------------------------------------
# # Filter out normal copy number segments
# # ------------------------------------------------------------------------------
amplificationCutoff = 30
gg = gg[cn>amplificationCutoff,]

#gg$set(y.field = 'cn')
gg$set(gr.labelfield = 'node.id')
gg$edges$mark(h=1) # This is how you change the 'curviness'
gg$edges$mark(lwd=log10(gg$edges$dt$AD_TUMOR))
plot(gg$gt, '2:10000000-100000000')

# # ------------------------------------------------------------------------------
# # Filter out reference edges if ALT edge is of high allele depth
# # ------------------------------------------------------------------------------

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
plot(gg$gt, '2:15000000-17000000')
plot(gg$gt, '2:10000000-100000000')

# # ------------------------------------------------------------------------------
# # Create strongly connected clusters
# # ------------------------------------------------------------------------------
# 
# # Get strong clusters
# gg$clusters(mode="strong")
# 
# # MYCN containing node
# MYCN_start = 16080683
# MYCN_end = 16087129
# mycn_node = gg$nodes[seqnames == "2" & start<= MYCN_start & end>=MYCN_start,]$dt$snode.id
# 
# # MYCN containing cluster
# mycn_cluster = gg$nodes[seqnames == "2" & start<= MYCN_start & end>=MYCN_start,]$dt$cluster
# 
# # MYCN containing strong cluster
# gg.sub = gg[cluster==mycn_cluster,]
# #View(gg.sub$dt)
# 
# # Filter Graph by Position of Nodes and edges
# #gg.sub = gg[seqnames == "2" & start > 10000000 & end < 100000000,] # before comma = nodes, after comma = edges
# 
# # Plot
# gg.sub$set(gr.labelfield = 'node.id')
# gg.sub$edges$mark(h=1) # This is how you change the 'curviness'
# gg.sub$edges$mark(lwd=log10(gg.sub$edges$dt$AD_CELLLINE))
# plot(gg.sub$gt, '2:10000000-100000000')
# 
# # Get node ID that contains MYCN (chr2:16,080,683-16,087,129 GRCh37/hg19)
# MYCN_start = 16080683
# MYCN_end = 16087129
# mycn_gg.sub_node = gg.sub$nodes[seqnames == "2" & start<= MYCN_start & end>=MYCN_start,]$dt$snode.id
# 
# # ------------------------------------------------------------------------------
# # Create walks object for neighborhood of MYCN
# # ------------------------------------------------------------------------------
# 
# #View(gg.sub$nodes$dt)
# #View(gg.sub$edges$dt)
# 
# # Create walks object
# w = gg.sub$walks() # takes 15 min for 22 nodes and 29 edges
# #save.image("/Volumes/Elements/nb-cl-wgs/playground/IMR575_withWalks.Rdata")
# #load("/Volumes/Elements/nb-cl-wgs/playground/IMR575_withWalks.Rdata")
# 
# 
# # subaset of circular walks that contain MYCN
# w.mycn = w[sapply(w$snode.id, function (x) (mycn_gg.sub_node %in% x) | (-mycn_gg.sub_node %in% x))] # only walks that contain this_node
# w.mycn = w.mycn[circular==TRUE]
# 
# # now, the very deep breakpoints need to be included
# deep_edge_threshold = 90
# 
# # TODO: add quality threshold etc.
# 
# # mark all edges larger than threshold
# gg.sub$edges$mark(deepEdge = gg.sub$edges$dt$AD_CELLLINE>=deep_edge_threshold)
# required_edges = data.frame("edgeID" = gg.sub$edges$dt$edge.id, "deepEdge" = gg.sub$edges$dt$deepEdge) %>% filter(deepEdge == TRUE) %>% .$edgeID
# idx = c()
# contains_req_edges = rep(FALSE,length(w.mycn))
# for (i in 1:length(w.mycn)){
#   this_walk_edges = abs(w.mycn$sedge.id[[i]])
#   if (setequal(required_edges, intersect(this_walk_edges, required_edges))) contains_req_edges[i] = TRUE
# }
# w.mycn$set(contains_req_edges = contains_req_edges)
# 
# w.mycn[contains_req_edges == TRUE]$dts() %>% View
# w.mycn.req_edges = w.mycn[contains_req_edges == TRUE]
# shortest_walk = w.mycn.req_edges[(order(wid))][1]
# 
# # Remove nodes until there is a shortes path that connects all exactly once.
# tes = w.mycn$snode.id[[1]]
# 
# hasdups = sapply(w.mycn$snode.id, function (sn) sum(duplicated(abs(sn)))>0)
# 
# w.mycn.nodups = w.mycn[sapply(w.mycn$snode.id, function (sn) sum(duplicated(abs(sn)))==0)]
# longest_walk_nodups = w.mycn.nodups[(order(wid, decreasing = T))][1]
# 

gg.sub = gg

# !sum(duplicated(abs(tes)))
# ------------------------------------------------------------------------------
# Proximity
# ------------------------------------------------------------------------------

# the interesting enhancer (from lowMYCN_aggregate_peaks): chr2 10532133-10534695
enhancer_of_interest = lowMYCN_aggregate_peaks %>%
  filter(seqnames == "chr2", start == 10532133)
seqlevelsStyle(enhancer_of_interest) = 'NCBI'
gff = readRDS(gzcon(url('http://mskilab.com/gGnome/hg19/gencode.v19.annotation.gtf.gr.rds')))
genes = gff %Q% (type == 'gene' & gene_name == "MYCN")
px = proximity(gg.sub, enhancer_of_interest, genes[, 'gene_name'], mc.cores=cores)
print(px)

# ------------------------------------------------------------------------------
# Overview Plot
# ------------------------------------------------------------------------------
library(BSgenome.Hsapiens.UCSC.hg19)

window = GRanges(
  seqnames = rep("2", 4),
  ranges = IRanges(start = c(10500000, 14600000, 32350000, 55100000),
                   end   = c(10700000, 16200000, 32400000, 55200000))
)

windows_with_chr = window
seqlevelsStyle(windows_with_chr) <- "UCSC" 

gTrackUtils.standardchrs = standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)[1:23]
gTrackUtils.binsize = 1000 
gTrackUtils.window_bins = tileGenome(seqlengths = seqlengths(BSgenome.Hsapiens.UCSC.hg19), tilewidth = gTrackUtils.binsize, cut.last.tile.in.chrom = T) %>%
  filter_by_overlaps(windows_with_chr, minoverlap = gTrackUtils.binsize/10)
seqlevels(gTrackUtils.window_bins, pruning.mode="coarse") = gTrackUtils.standardchrs
make_bigwig_gTrack = function(fname, name = "NoName", y0=0, y1=50){
  this_track = read_bigwig(fname, overlap_ranges = windows_with_chr)
  seqlevels(this_track, pruning.mode="coarse") = gTrackUtils.standardchrs
  this_track = binnedAverage(gTrackUtils.window_bins, coverage(this_track, weight=this_track$score), "score")
  this_gt = gTrack(this_track, y.field="score", bar=TRUE, name = name, y0=y0, y1=y1) #, 
  return(this_gt)
}

Sys.setenv(GENCODE_DIR = "/Volumes/Elements/")
gff = readRDS(gzcon(url('http://mskilab.com/gGnome/hg19/gencode.v19.annotation.gtf.gr.rds')))
genes = gff %Q% (type == 'gene' & gene_name == "MYCN")
gencode = track.gencode(stack.gap = 2e5, cex.label = 0.5, height = 20, name = 'GENCODE', build='hg19', labels.suppress.gr=TRUE, gene.collapse = TRUE, genes = gff %Q% (type == 'gene' & gene_type == "protein_coding") %>% .$gene_name)

gt.copynumber = gTrack(ascat_cnv_gr %>% filter(Sample == "NBL26"), y.field = "TumorTotalCopyNumber")
gt.NoMYCNExpr_H3K27ac = make_bigwig_gTrack("/Volumes/Elements/MYCNAmplicon/Results/Boeva_nonMYCNExpr_H3K27ac_MeanFC.bw", name = "Non MYCN Expr NB (H3K27ac FC)", y0=0, y1=10)
gt.LowMYCNExpr_H3K27ac = make_bigwig_gTrack("/Volumes/Elements/MYCNAmplicon/Results/Boeva_lowMYCNExpr_H3K27ac_MeanFC.bw", name = "Low MYCN Expr NB (H3K27ac FC)", y0=0, y1=10)
gt.CLBGA_GATA3 = make_bigwig_gTrack("/Volumes/Elements/nb-cl-chipseq-results/bam/Boeva_CLB-GA_GATA3.trimmed.bwa_hg19.rmdup.filtered.bw", name = "CLB-GA GATA3", y0=0, y1=5)
gt.CLBGA_PHOX2B = make_bigwig_gTrack("/Volumes/Elements/nb-cl-chipseq-results/bam/Boeva_CLB-GA_PHOX2B.trimmed.bwa_hg19.rmdup.filtered.bw", name = "CLB-GA PHOX2B", y0=0, y1=5)
gt.CLBGA_HAND2 = make_bigwig_gTrack("/Volumes/Elements/nb-cl-chipseq-results/bam/Boeva_CLB-GA_HAND2.trimmed.bwa_hg19.rmdup.filtered.bw", name = "CLB-GA HAND2", y0=0, y1=5)

noMYCN_aggregate_peaks = read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_noMYCNExpr_H3K27ac_AggregatePeaks.bed")
noMYCN_aggregate_peaks.gt = gTrack(noMYCN_aggregate_peaks, height = 10)

lowMYCN_aggregate_peaks = read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_lowMYCNExpr_H3K27ac_AggregatePeaks.bed")
lowMYCN_aggregate_peaks.gt = gTrack(lowMYCN_aggregate_peaks, height = 10)

# Create gTrack for all non-blacklisted rearrangements
svaba_forplotting = jJ(rafile=out_vcf, geno = TRUE) # I needed to change utils.R in the gGnome packages (VariantAnnotation::geno instead of geno)
blacklisted_svaba_forplotting_indices = sort(unique(c(
  queryHits(findOverlaps(svaba_forplotting$left, blacklist_gr)), 
  queryHits(findOverlaps(svaba_forplotting$right, blacklist_gr))
)))
svaba_forplotting = svaba_forplotting[-blacklisted_svaba_forplotting_indices]

gg_forplotting = gG(juncs=svaba_forplotting)
gg_forplotting$edges$mark(h=3) # This is how you change the 'curviness'
gg_forplotting$edges$mark(lwd=0.0001*(gg_forplotting$edges$dt$AD_CELLLINE))
sv.gt = gg_forplotting$gt
sv.gt$height = 50
plot(sv.gt, window)

proximity.gt = px$gt
proximity.gt$height = 10

pdf("/Volumes/Elements/MYCNAmplicon/Results/NBL26_Overview.pdf", height=7, width=7, onefile=F, useDingbats = F)
plot(c(gencode, gt.copynumber, sv.gt, proximity.gt, gt.NoMYCNExpr_H3K27ac, noMYCN_aggregate_peaks.gt, gt.LowMYCNExpr_H3K27ac, lowMYCN_aggregate_peaks.gt, gt.CLBGA_PHOX2B, gt.CLBGA_GATA3, gt.CLBGA_HAND2), window,
     y.grid.col=F, y.grid.lwd=F, sep.lwd=0, sep.col=F, yaxis.pretty = 1)
dev.off()




