library(ggbio)
library(plyranges)
library(BSgenome.Hsapiens.UCSC.hg19)

window = GRanges(
  seqnames = rep("chr2", 3),
  ranges = IRanges(start = c(11500000, 15000000, 15500000),
                   end   = c(13000000, 15100000, 16500000))
)

ggbio(trackWidth = 10, buffer = 0, radius = 10) + 
  circle(window, geom = "ideo", fill = "gray70") +
  circle(window, geom = "scale", size = 2) +
  circle(window, geom = "text", aes(label = seqnames), vjust = 0, size = 3)

library(circlize)
tp_family = readRDS(system.file(package = "circlize", "extdata", "tp_family_df.rds"))
head(tp_family)
circos.genomicInitialize(tp_family, tickLabelsStartFromZero=F)
circos.track(ylim = c(0, 1), 
             bg.col = c("#FF000040", "#00FF0040", "#0000FF40"), 
             bg.border = NA, track.height = 0.05)

circos.clear()
#window_df = as.data.frame(window)
#window_df$seqnames = paste0(1:nrow(window_df), "_", window_df$seqnames)
circos.genomicInitialize(as.data.frame(window), tickLabelsStartFromZero = F)

hg19sub = keepStandardChromosomes(GRangesForBSGenome("hg19"), pruning.mode = "coarse")
bins <- tileGenome(seqinfo(hg19sub), tilewidth=1000000, cut.last.tile.in.chrom=TRUE)
mycols = rand_color(n=100, luminosity = "dark")
bgcolor = "gray98"

cnv = read.table("/Volumes/Elements/nb-cl-wgs/nb-cl-wgs-freec/CHP-212/CHP-212_IlluminaWGS.hg19.rmdup.fixRG.bam_CNVsfillgaps.txt",
                 header=F, sep="\t")
colnames(cnv) = c("seqnames", "start", "end", "width", "strand", "copy.number", "copy.number.call")
cnv = makeGRangesFromDataFrame(cnv, keep.extra.columns = T)
seqlevels(cnv,pruning.mode="coarse") = seqlevels(hg19sub)
seqlengths(cnv) = seqlengths(hg19sub)

cnv = binnedAverage(bins=bins, coverage(cnv, weight=cnv$copy.number), "y") %>% as.data.frame() %>% dplyr::select(seqnames, start, end, y)

circos.genomicTrack(as.data.frame(cnv), 
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value[[1]],  ytop.column = 1, ybottom = 0,
                                         col = ifelse(value[[1]] > 0, "firebrick4", ifelse(value[[1]] < 0, "steelblue4", "lightgray")),
                                         border = ifelse(value[[1]] > 0, "firebrick4", ifelse(value[[1]] < 0, "steelblue4", "lightgray")), ...)
                    },
                    track.height=0.5*circos.par("track.height"),
                    bg.border=NA, bg.col = bgcolor)


circos.initializeWithIdeogram(species = "hg19")
circos.genomicTrack(as.data.frame(cnv), 
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value[[1]],  ytop.column = 1, ybottom = 0,
                                         col = ifelse(value[[1]] > cnv_normal, "firebrick4", ifelse(value[[1]] < cnv_normal, "steelblue4", "lightgray")),
                                         border = ifelse(value[[1]] > cnv_normal, "firebrick4", ifelse(value[[1]] < cnv_normal, "steelblue4", "lightgray")), ...)
                    },
                    track.height=0.5*circos.par("track.height"),
                    bg.border=NA, bg.col = "white")

