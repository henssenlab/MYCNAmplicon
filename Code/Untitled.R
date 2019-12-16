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
