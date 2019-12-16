# library(gTrack)
# devtools::load_all("/Volumes/Elements/MYCNAmplicon/Ext/gGnome-master") # includes my changes
# library(BSgenome.Hsapiens.UCSC.hg19)
# library(ggplot2)
# library(dplyr)
# library(parallel)
# library(plyranges)
library(BSgenome.Hsapiens.UCSC.hg19)
gTrackUtils.standardchrs = standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)[1:23]
gTrackUtils.binsize = 1000 
gTrackUtils.window_bins = tileGenome(seqlengths = seqlengths(BSgenome.Hsapiens.UCSC.hg19),tilewidth = gTrackUtils.binsize, cut.last.tile.in.chrom = T) %>%
  filter_by_overlaps(window, minoverlap = gTrackUtils.binsize/10)
seqlevels(gTrackUtils.window_bins, pruning.mode="coarse") = gTrackUtils.standardchrs
make_bigwig_gTrack = function(fname, name = "NoName", y0=0, y1=50){
  this_track = read_bigwig(fname, overlap_ranges = window)
  seqlevels(this_track , pruning.mode="coarse") = gTrackUtils.standardchrs
  this_track = binnedAverage(gTrackUtils.window_bins, coverage(this_track, weight=this_track$score), "score")
  this_gt = gTrack(this_track, y.field="score", bar=TRUE, name = name, y0=y0, y1=y1) #, 
  return(this_gt)
}