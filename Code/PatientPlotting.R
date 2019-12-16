library(gTrack)
devtools::load_all("/Volumes/Elements/MYCNAmplicon/Ext/gGnome-master") # includes my changes
#devtools::load_all("/Volumes/Elements/nb-reconstruct-amplicon/scripts/gTrack-master") # includes my changes
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)
library(dplyr)
library(parallel)
library(plyranges)

load("/Volumes/Transcend/PalmTrees/Analysis/WorkspaceData/ascat_cnv.Rdata")

blacklist_fname = "~/Downloads/hg19-blacklist.v2.bed"
blacklist = read.table(blacklist_fname, sep="\t", header=F)
colnames(blacklist) = c("Chr", "Start", "End", "Class")
blacklist$Chr = gsub("chr", "", blacklist$Chr)
blacklist_gr = makeGRangesFromDataFrame(blacklist, keep.extra.columns = T)

Sys.setenv(GENCODE_DIR = "/Volumes/Elements/")
gff = readRDS(gzcon(url('http://mskilab.com/gGnome/hg19/gencode.v19.annotation.gtf.gr.rds')))
genes = gff %Q% (type == 'gene' & gene_name == "MYCN")
seqlevels(genes) = paste0("chr", seqlevels(genes))
gencode = track.gencode(stack.gap = 2e5, cex.label = 0.5, height = 15, name = 'GENCODE', build='hg19', labels.suppress.gr=TRUE, gene.collapse = TRUE, genes = gff %Q% (type == 'gene' & gene_type == "protein_coding") %>% .$gene_name)

se = read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_lowMYCN_CRCdriven_SE.bed")
strand(se) = "*"
gt.se = gTrack(se, col="gold2", height = 5)

mycn_enh = read_bed("/Volumes/Elements/MYCNAmplicon/Results/Boeva_nMNA_MYCN_Enhancers.bed")
strand(mycn_enh) = "*"
gt.mycn_enh = gTrack(mycn_enh, col="firebrick3", height = 5)


gt.aggregateh3k27ac = gTrack("/Volumes/Elements/MYCNAmplicon/Results/Boeva_lowMYCNExpr_H3K27ac_MeanFC.bw", bar=TRUE, name = "Aggr. H3K27ac", y0=0, y1=15)

# ------------------------------------------------------------------------------
# NB2001
# ------------------------------------------------------------------------------

window = GRanges(
  seqnames = c("1", "1", "2", "2", "2", "2", "2"),
  ranges = IRanges(start = c(30500000, 32200000, 2350000, 6250000, 15500000, 36500000, 58000000),
                   end   = c(31000000, 32500000, 2550000, 7000000, 16500000, 37000000, 58500000))
)

out_vcf = "/Volumes/Elements/MYCNAmplicon/Data/SvABA/NB2001/NB2001.svaba.unfiltered.somatic.sv.vcf.chr1_chr2.vcf" # needs to be a vcf
svaba = jJ(rafile=out_vcf, geno = TRUE) # I need to change utils.R in the gGnome packages (VariantAnnotation::geno instead of geno)
blacklisted_svaba_indices = sort(unique(c(
  queryHits(findOverlaps(svaba$left, blacklist_gr)), # those are the filtered bins
  queryHits(findOverlaps(svaba$right, blacklist_gr))
)))
svaba_filtered = svaba[-blacklisted_svaba_indices]
svaba_filtered = svaba_filtered[FILTER == "PASS"]
gg = gG(juncs=svaba_filtered)
gg$edges$mark(lwd=0.1*(gg$edges$dt$AD_TUMOR))
gt.sv = gg$gt
gt.sv$height = 50
gt.copynumber = gTrack(ascat_cnv_gr %>% filter(Sample == "NB2001"), y.field = "TumorTotalCopyNumber")

pdf("/Volumes/Elements/MYCNAmplicon/Results/SVPlot_NB2001.pdf", height=7, width=7, useDingbats = F, onefile=F)
plot(c(gencode, gt.se, gt.mycn_enh, gt.aggregateh3k27ac, gt.copynumber, gt.sv), window)
dev.off()

# ------------------------------------------------------------------------------
# NBL26
# ------------------------------------------------------------------------------

window = GRanges(
  seqnames = rep("2", 4),
  ranges = IRanges(start = c(10500000, 14600000, 32350000, 55100000),
                   end   = c(10700000, 16500000, 32400000, 55200000))
)

out_vcf = "/Volumes/Elements/MYCNAmplicon/Data/SvABA/NBL26/NBL26.svaba.unfiltered.somatic.sv.vcf.chr2_10M_75M.vcf" # needs to be a vcf
svaba = jJ(rafile=out_vcf, geno = TRUE) # I need to change utils.R in the gGnome packages (VariantAnnotation::geno instead of geno)
blacklisted_svaba_indices = sort(unique(c(
  queryHits(findOverlaps(svaba$left, blacklist_gr)), # those are the filtered bins
  queryHits(findOverlaps(svaba$right, blacklist_gr))
)))
svaba_filtered = svaba[-blacklisted_svaba_indices]
svaba_filtered = svaba_filtered[FILTER == "PASS"]
gg = gG(juncs=svaba_filtered)
gg$edges$mark(lwd=0.1*(gg$edges$dt$AD_TUMOR))
gt.sv = gg$gt
gt.sv$height = 50
gt.copynumber = gTrack(ascat_cnv_gr %>% filter(Sample == "NBL26"), y.field = "TumorTotalCopyNumber")
pdf("/Volumes/Elements/MYCNAmplicon/Results/SVPlot_NBL26.pdf", height=7, width=7, useDingbats = F, onefile=F)
plot(c(gencode, gt.se, gt.mycn_enh, gt.aggregateh3k27ac, gt.copynumber, gt.sv), window)
dev.off()

# ------------------------------------------------------------------------------
# NBL32
# ------------------------------------------------------------------------------

window = GRanges(
  seqnames = rep("2", 4),
  ranges = IRanges(start = c(15100000, 16800000, 17400000, 18000000),
                   end   = c(16500000, 17100000, 17800000, 18250000))
)

out_vcf = "/Volumes/Elements/MYCNAmplicon/Data/SvABA/NBL32/NBL32.svaba.unfiltered.somatic.sv.vcf.chr2_10M_75M.vcf" # needs to be a vcf
gt.copynumber = gTrack(ascat_cnv_gr %>% filter(Sample == "NBL32"), y.field = "TumorTotalCopyNumber")

svaba = jJ(rafile=out_vcf, geno = TRUE) # I need to change utils.R in the gGnome packages (VariantAnnotation::geno instead of geno)
blacklisted_svaba_indices = sort(unique(c(
  queryHits(findOverlaps(svaba$left, blacklist_gr)), # those are the filtered bins
  queryHits(findOverlaps(svaba$right, blacklist_gr))
)))
svaba_filtered = svaba[-blacklisted_svaba_indices]
svaba_filtered = svaba_filtered[FILTER == "PASS"]
gg = gG(juncs=svaba_filtered)
gg$edges$mark(lwd=0.025*(gg$edges$dt$AD_TUMOR))
gt.sv = gg$gt
gt.sv$height = 50
pdf("/Volumes/Elements/MYCNAmplicon/Results/SVPlot_NBL32.pdf", height=7, width=7, useDingbats = F, onefile=F)
plot(c(gencode, gt.se, gt.mycn_enh, gt.aggregateh3k27ac, gt.copynumber, gt.sv), window)
dev.off()

