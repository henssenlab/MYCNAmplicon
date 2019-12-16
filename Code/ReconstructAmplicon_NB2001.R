library(gTrack)
devtools::load_all("/Volumes/Elements/MYCNAmplicon/Ext/gGnome-master") # includes my changes
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
wgs_cov_cov = ascat_cnv_gr %>% filter(Sample == "NB2001") %>%
  coverage(weight = "TumorTotalCopyNumber")

window = GRanges(
  seqnames = c("1", "2", "2", "2", "2", "2"),
  ranges = IRanges(start = c(30500000, 2350000, 6000000, 15500000, 36500000, 58000000),
                   end   = c(32500000, 2550000, 7000000, 16500000, 37500000, 58500000))
)

# TODO filter bw files
#wgs_cov = plyranges::read_bigwig("/Volumes/Elements/nb-cl-wgs/nb-cl-wgs-bam/IMR-5-75_IlluminaWGS.hg19.rmdup.fixRG.bw")
#wgs_cov_cov = coverage(wgs_cov, weight=wgs_cov$score)

# ------------------------------------------------------------------------------
# Prepare svaba vcf
# ------------------------------------------------------------------------------

# Filter vcf files 
# changed "ID=PL,Number=." to "ID=PL,Number=G" in header
# changed column names to tumor and normal
# bgzip somatic_run.svaba.unfiltered.somatic.sv.vcf
# bcftools index somatic_run.svaba.unfiltered.somatic.sv.vcf.gz
# bcftools view somatic_run.svaba.unfiltered.somatic.sv.vcf.gz --regions 1,2 > NB2001.svaba.unfiltered.somatic.sv.vcf.chr1_chr2.vcf

# ------------------------------------------------------------------------------
# Create gGraph from junctions 
# ------------------------------------------------------------------------------

out_vcf = "/Volumes/Elements/MYCNAmplicon/Data/SvABA/NB2001/NB2001.svaba.unfiltered.somatic.sv.vcf.chr1_chr2.vcf" # needs to be a vcf

svaba = jJ(rafile=out_vcf, geno = TRUE) # I need to change utils.R in the gGnome packages (VariantAnnotation::geno instead of geno)
blacklisted_svaba_indices = sort(unique(c(
  queryHits(findOverlaps(svaba$left, blacklist_gr)), # those are the filtered bins
  queryHits(findOverlaps(svaba$right, blacklist_gr))
)))
svaba_filtered = svaba[-blacklisted_svaba_indices]
svaba_filtered = svaba_filtered[AD_TUMOR>10] 

gg = gG(juncs=svaba_filtered)

#gg = gg[seqnames == "chr2"]
#gg = gg[start<33141000 | start>33142000,]

# Assign each node with its mean copy number signal
load("/Volumes/Transcend/PalmTrees/Analysis/WorkspaceData/ascat_cnv.Rdata")
seqlevelsStyle(ascat_cnv_gr) <- "NCBI" # get rid of chr in seqlevels
seqlevels(ascat_cnv_gr, pruning.mode="coarse") = seqlevels(gg$nodes$gr)
seqlengths(ascat_cnv_gr) = seqlengths(gg$nodes$gr)
wgs_cov_cov = ascat_cnv_gr %>% filter(Sample == "NB2001") %>%
  coverage(weight = "TumorTotalCopyNumber")
gg$nodes$mark(cn=binnedAverage(gg$nodes$gr, wgs_cov_cov, "cn", na.rm=T) %>% .$cn)
gg$nodes$mark(cn_log10=log10(gg$nodes$dt$cn + 0.0001))

#lapply(seq_along(gg$nodes$gr), function(i) wgs_cov %>% filter_by_overlaps(gg$nodes$gr[i]) %>% .$score %>% median)

#gg$nodes$dt %>% View
#gg$set(y.field = 'cn_log10')
gg$set(gr.labelfield = 'node.id')
gg$edges$mark(h=5) # This is how you change the 'curviness'
#gg$edges$mark(lwd=0.1*(gg$edges$dt$AD_TUMOR))
plot(gg$gt, window)

# ------------------------------------------------------------------------------
# Filter out normal copy number segments
# ------------------------------------------------------------------------------
amplificationCutoff = 30
gg = gg[cn>amplificationCutoff,]

#gg$set(y.field = 'cn')
gg$set(gr.labelfield = 'node.id')
gg$edges$mark(h=1) # This is how you change the 'curviness'
gg$edges$mark(lwd=log10(gg$edges$dt$AD_TUMOR))
plot(gg$gt, '2:1-100000')

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
gg$edges$mark(lwd=log10(gg$edges$dt$AD_TUMOR))
plot(gg$gt, '2:14000000-20000000')

# ------------------------------------------------------------------------------
# Create strongly connected clusters
# ------------------------------------------------------------------------------

# Get strong clusters
gg$clusters(mode="strong")

# MYCN containing node
MYCN_start = 16080683
MYCN_end = 16087129
mycn_node = gg$nodes[seqnames == "2" & start<= MYCN_start & end>=MYCN_start,]$dt$snode.id

# MYCN containing cluster
mycn_cluster = gg$nodes[seqnames == "2" & start<= MYCN_start & end>=MYCN_start,]$dt$cluster

# MYCN containing strong cluster
gg.sub = gg[cluster==mycn_cluster,]
#View(gg.sub$dt)

# Filter Graph by Position of Nodes and edges
#gg.sub = gg[seqnames == "2" & start > 10000000 & end < 100000000,] # before comma = nodes, after comma = edges

# Plot
gg.sub$set(gr.labelfield = 'node.id')
gg.sub$edges$mark(h=1) # This is how you change the 'curviness'
gg.sub$edges$mark(lwd=log10(gg.sub$edges$dt$AD_CELLLINE))
plot(gg.sub$gt, '2:14000000-20000000')

# Get node ID that contains MYCN (chr2:16,080,683-16,087,129 GRCh37/hg19)
MYCN_start = 16080683
MYCN_end = 16087129
mycn_gg.sub_node = gg.sub$nodes[seqnames == "2" & start<= MYCN_start & end>=MYCN_start,]$dt$snode.id

# ------------------------------------------------------------------------------
# Create walks object for neighborhood of MYCN
# ------------------------------------------------------------------------------

#View(gg.sub$nodes$dt)
#View(gg.sub$edges$dt)

# Create walks object
w = gg.sub$walks() # takes not long
#save.image("/Volumes/Elements/MYCNAmplicon/Figures/NBL32Reconstruction_PLAYGROUND.Rdata")
#load("/Volumes/Elements/MYCNAmplicon/Figures/NBL32Reconstruction_PLAYGROUND.Rdata")


# subset of circular walks that contain MYCN
w.mycn = w[sapply(w$snode.id, function (x) (mycn_gg.sub_node %in% x) | (-mycn_gg.sub_node %in% x))] # only walks that contain this_node
w.mycn = w.mycn[circular==TRUE]

# now, the very deep breakpoints need to be included
deep_edge_threshold = 90

# TODO: add quality threshold etc.

# mark all edges larger than threshold
gg.sub$edges$mark(deepEdge = gg.sub$edges$dt$AD_CELLLINE>=deep_edge_threshold)
required_edges = data.frame("edgeID" = gg.sub$edges$dt$edge.id, "deepEdge" = gg.sub$edges$dt$deepEdge) %>% filter(deepEdge == TRUE) %>% .$edgeID
idx = c()
contains_req_edges = rep(FALSE,length(w.mycn))
for (i in 1:length(w.mycn)){
  this_walk_edges = abs(w.mycn$sedge.id[[i]])
  if (setequal(required_edges, intersect(this_walk_edges, required_edges))) contains_req_edges[i] = TRUE
}
w.mycn$set(contains_req_edges = contains_req_edges)
w.mycn[contains_req_edges == TRUE]$dts() %>% View
w.mycn.req_edges = w.mycn[contains_req_edges == TRUE]
shortest_walk = w.mycn.req_edges[(order(wid))][1]

# Remove nodes until there is a shortes path that connects all exactly once.
tes = w.mycn$snode.id[[1]]
hasdups = sapply(w.mycn$snode.id, function (sn) sum(duplicated(abs(sn)))>0)

w.mycn = w.mycn[w.mycn$circular]
w.mycn.nodups = w.mycn[sapply(w.mycn$snode.id, function (sn) sum(duplicated(abs(sn)))==0)]
longest_walk_nodups = w.mycn.nodups[(order(wid, decreasing = T))][1]

se = read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_lowMYCN_CRCdriven_SE.bed")
strand(se) = "*"
gt.se = gTrack(se, col="gold2", height = 5)

!sum(duplicated(abs(tes)))
# ------------------------------------------------------------------------------
# Proximity analysis
# ------------------------------------------------------------------------------

enhancers = plyranges::read_bed("/Volumes/Elements/MYCNAmplicon/Figures/Boeva_nMNA_MYCN_Enhancers.bed")
seqlevelsStyle(enhancers) <- "NCBI" # get rid of chr in seqlevels
enhancers = enhancers %>% filter(name == "e4")
gff = readRDS(gzcon(url('http://mskilab.com/gGnome/hg19/gencode.v19.annotation.gtf.gr.rds')))
genes = gff %Q% (type == 'gene' & gene_name == "MYCN")
px = proximity(gg.sub, enhancers, genes[, 'gene_name'], mc.cores=cores)

# ------------------------------------------------------------------------------
# Overview Plot
# ------------------------------------------------------------------------------

window = GRanges(
  seqnames = c("1", "2", "2", "2", "2", "2"),
  ranges = IRanges(start = c(30500000, 2350000, 6000000, 15500000, 36500000, 58000000),
                   end   = c(32500000, 2550000, 7000000, 16500000, 37500000, 58500000))
)

gff = readRDS(gzcon(url('http://mskilab.com/gGnome/hg19/gencode.v19.annotation.gtf.gr.rds')))
genes = gff %Q% (type == 'gene' & gene_name == "MYCN")
Sys.setenv(GENCODE_DIR = "/Volumes/Elements/")
gencode = track.gencode(stack.gap = 2e5, cex.label = 0.5, height = 20, name = 'GENCODE', build='hg19', labels.suppress.gr=TRUE, gene.collapse = TRUE, genes = gff %Q% (type == 'gene' & gene_type == "protein_coding") %>% .$gene_name)

gt.copynumber = gTrack(ascat_cnv_gr %>% filter(Sample == "NB2001"), y.field = "TumorTotalCopyNumber")
gt.NoMYCNExpr_H3K27ac = gTrack("/Volumes/Elements/MYCNAmplicon/Results/Boeva_nonMYCNExpr_H3K27ac_MeanFC.bw", bar=TRUE, name = "Non MYCN Expr NB (H3K27ac FC)", y0=0, y1=20)
#plot(gt.NoMYCNExpr_H3K27ac, window, y.grid.col=F, sep.lwd=F)
gt.LowMYCNExpr_H3K27ac = gTrack("/Volumes/Elements/MYCNAmplicon/Results/Boeva_lowMYCNExpr_H3K27ac_MeanFC.bw", bar=TRUE, name = "Low MYCN Expr NB (H3K27ac FC)", y0=0, y1=20)
gt.CLBGA_GATA3 = gTrack("/Volumes/Elements/nb-cl-chipseq-results/bam/Boeva_CLB-GA_GATA3.trimmed.bwa_hg19.rmdup.filtered.bw", bar=TRUE, name = "CLB-GA GATA3", y0=0, y1=10)
gt.CLBGA_PHOX2B = gTrack("/Volumes/Elements/nb-cl-chipseq-results/bam/Boeva_CLB-GA_PHOX2B.trimmed.bwa_hg19.rmdup.filtered.bw", bar=TRUE, name = "CLB-GA PHOX2B", y0=0, y1=10)
gt.CLBGA_HAND2 = gTrack("/Volumes/Elements/nb-cl-chipseq-results/bam/Boeva_CLB-GA_HAND2.trimmed.bwa_hg19.rmdup.filtered.bw", bar=TRUE, name = "CLB-GA HAND2", y0=0, y1=10)

gg$set(gr.labelfield = F)

pdf("/Volumes/Elements/MYCNAmplicon/Results/NBL2001_Overview.pdf", height=10, width=10, onefile=F, useDingbats = F)
plot(c(gencode, gg$gt, gt.copynumber, gt.LowMYCNExpr_H3K27ac, gt.NoMYCNExpr_H3K27ac, gt.CLBGA_PHOX2B, gt.CLBGA_GATA3, gt.CLBGA_HAND2), window, y.grid.col=F, y.grid.lwd=F, sep.lwd=0, sep.col=F)
dev.off()
