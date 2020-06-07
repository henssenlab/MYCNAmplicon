rm(list=ls())
library(Rcpp)
library(dplyr)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg19)
library(RcppCNPy)
library(data.table)
library(rtracklayer)
sourceCpp("/Volumes/Transcend/HiC/straw-master/R/straw-R.cpp")

# MYCN chr2:16,080,683-16,087,129
# MYCN TSS chr2:16,080,683
# MYCN TSS BIN 5kb 2:16080000-2:16085000 or 2:16075000-2:16085000

virtual4C = function(hic_fname, viewpoint_chr, viewpoint_start, viewpoint_end, bw_fname, binsize=5000, normalization_type = "KR"){

  # hic_fname = "/Volumes/Transcend/HiC/Normalization/Raw/IMR575.hic"
  # hic_fname = "/Volumes/Transcend/HiC/Normalization/Raw/NGP.hic"
  # hic_fname = "/Volumes/Elements/MariaSalaDaten/HiC/IMR5/IMR5_hg19_canonical_MAPQ30_merged.hic"
  # bw_fname = "~/Desktop/test.bw"
  # viewpoint_chr = "chr2"
  # viewpoint_chr_index = which(viewpoint_chr == seqnames(BSgenome.Hsapiens.UCSC.hg19))
  # viewpoint_start = 16075000 # this is the actual start; only reads/info after this coordinate
  # viewpoint_end = 16085000 # this is the actual end; only reads/info before this coordinate
  # normalization_type = "KR"
  # binsize = 5000
  
  viewpoint_all_bins = seq(viewpoint_start, viewpoint_end, binsize)
  #chromosomes = standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)
  #chr_lengths = seqlengths(BSgenome.Hsapiens.UCSC.hg19)[standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)]
  #which("chr2" == seqnames(BSgenome.Hsapiens.UCSC.hg19))
  
  chromosomes = viewpoint_chr
  chr_lengths = seqlengths(BSgenome.Hsapiens.UCSC.hg19)[viewpoint_chr]
  
  contacts = list()
  for (this_chr in chromosomes){
  
    empty_template = 
      data.frame(
        "Pos" = rep(seq(0, chr_lengths[[this_chr]], binsize), length(viewpoint_all_bins)),
        "ViewpointPos" = rep(viewpoint_all_bins, each=ceiling(chr_lengths[[this_chr]] / binsize)),
        "value" = 0)

    if (viewpoint_chr == this_chr){
      this_pair = straw_R(paste0(normalization_type, " ", hic_fname, " ", gsub("chr", "", this_chr), " ", gsub("chr", "", viewpoint_chr), ":", sprintf("%.0f", viewpoint_start), ":", sprintf("%.0f", viewpoint_end), " BP ", sprintf("%.0f", binsize))) %>% 
        as_tibble() %>%
        filter(!is.na(counts))
      this_pair_rev = this_pair # copy the data
      colnames(this_pair_rev) = c("y", "x", "counts") # flip x and y to get all values in x and all values in y, easy to code but honestly too memory intensive
      this_pair_db = bind_rows(this_pair, this_pair_rev)
      this_pair_db = this_pair_db %>% dplyr::rename(Pos = x, ViewpointPos = y, value = counts)
      this_pair_db = this_pair_db %>% filter(ViewpointPos >= viewpoint_start, ViewpointPos <= viewpoint_end)
    } else {
      stop("Interchromosomal data not implemented yet")
    }
    
    contacts[[this_chr]] = 
      this_pair_db %>%
      bind_rows(empty_template) %>%
      group_by(Pos, ViewpointPos) %>%
      summarise(value=sum(value)) %>%
      ungroup() %>%
      group_by(Pos) %>%
      summarise(value = mean(value)) %>%
      ungroup() %>%
      mutate(Chr = this_chr) %>%
      dplyr::select(Chr, Pos, value)
  }
  contacts = do.call(rbind, contacts)
  
  contacts_gr = contacts %>% 
    mutate(Start=Pos, End = Pos+(binsize-1), Name = ".", score=value) %>%
    dplyr::select(Chr, Start, End, score) %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)
  seqlengths(contacts_gr) = seqlengths(BSgenome.Hsapiens.UCSC.hg19)[seqlevels(contacts_gr)]
  rtracklayer::export.bw(contacts_gr, bw_fname)
  
}

hic_fnames = c(
  "/Volumes/Elements/MariaSalaDaten/HiC/CHP/CHP_hg19_canonical_MAPQ30_merged.hic",
  "/Volumes/Elements/MariaSalaDaten/HiC/IMR5/IMR5_hg19_canonical_MAPQ30_merged.hic",
  "/Volumes/Elements/MariaSalaDaten/HiC/KELLY-mearged/KELLY_hg19_canonical_MAPQ30_merged.hic",
  "/Volumes/Elements/MariaSalaDaten/HiC/Lan/Lan_hg19_canonical_MAPQ30_merged.hic",
  "/Volumes/Transcend/HiC/NGP.hic", 
  "/Volumes/Transcend/HiC/SKNDZ/4DNFIL1FQDXE.hic")
  
for (hic_fname in hic_fnames){
  virtual4C(
    hic_fname = hic_fname, 
    viewpoint_chr = "chr2", viewpoint_start=16075000, viewpoint_end=16085000, 
    bw_fname = paste0(hic_fname, ".MYCN4C_Aug28.Normalization_KR_Viewpoint_16075000_16085000_Res_5000.bw"), 
    binsize=5000, normalization_type = "KR")
}
