rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(plyranges)
library(data.table)
library(biomaRt)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggrepel)
library(ggforce)
source("/Volumes/Elements/MYCNAmplicon/Code/hasOverlap.R")
source("/Volumes/Elements/MYCNAmplicon/Code/CustomThemes.R")

meta = fread("/Volumes/Elements/MYCNAmplicon/Data/TCGA_ASCAT/summary.ascatTCGA.penalty70.txt") %>% 
  as_tibble() %>%
  filter(pass == TRUE, rep == TRUE) # pass quality control

nrow(meta)
length(unique(meta$name))

seg = list()
for (i in 1:nrow(meta)){
  seg[[i]] = 
    fread(paste0("/Volumes/Elements/MYCNAmplicon/Data/TCGA_ASCAT/segments/", meta[[i, "name"]], ".segments.txt")) %>%
    as_tibble() 
}
seg = do.call(rbind, seg)

cases_per_entity = 
  seg %>%  
  inner_join(meta %>% dplyr::select(name, cancer_type) %>% dplyr::rename(sample = name)) %>% 
  group_by(cancer_type) %>%
  summarise(n = n_distinct(sample))

seg = seg %>% 
  filter(nMajor >= 9) %>%
  inner_join(meta %>% dplyr::select(name, cancer_type) %>% dplyr::rename(sample = name))

# ------------------------------------------------------------------------------
# how many MYC(N) amplifizierte Fälle?
# ------------------------------------------------------------------------------

# 5         MYC               8      128747680    128753674
# 6        MYCN               2       16080686     16087129


seg %>%
  mutate(MYCamp =  hasOverlap_withChr(chr, startpos, endpos, 8, 128747680, 128753674),
         MYCNamp = hasOverlap_withChr(chr, startpos, endpos, 2, 16080686, 16087129)) %>%
  group_by(sample) %>%
  summarise(MYCamp = any(MYCamp),
            MYCNamp = any(MYCNamp),
            cancer_type = cancer_type[1]) %>%
  ungroup() %>% 
  group_by(cancer_type) %>%
  summarise(
            MYCamp_n = sum(MYCamp),
            MYCNamp_n = sum(MYCNamp),
          ) %>%
  ungroup() %>%
  full_join(cases_per_entity) %>%
  mutate(MYCamp_perc = 100 * MYCamp_n / n,
         MYCNamp_perc = 100 * MYCNamp_n / n) %>% 
  print()
  
entities = unique(seg$cancer_type)
entities = c("ESCA", "BRCA", "LUSC", "GBM", "SARC", "LGG")
entities = c("ESCA")

#hgnc_symbols = c("MYCN", "MYC", "EGFR", "CDK4", "CCND1", "TERC", "TERT", "MDM2", "FRS2", "HMGA2", "MCL1")
#hgnc_symbols = c("CDK4", "MDM2", "FRS2", "HMGA2")
#hgnc_symbols = c("PDGFRA")
#hgnc_symbols = c("MYC", "MYCN", "CDK4", "MDM2", "CCND1", "EGFR", "PDGFRA")
hgnc_symbols = c("MYC")
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl", host="grch37.ensembl.org")
genes =  getBM(attributes=c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
               filters="hgnc_symbol", 
               values=hgnc_symbols, 
               mart=ensembl)

for (gene_idx in 1:nrow(genes)){
  
  # Get list of samples that harbor amplification of the gene of interest
  ampl_samples = seg %>%
    filter(hasOverlap_withChr(chr, startpos, endpos, 
                              genes[gene_idx, "chromosome_name"],
                              genes[gene_idx, "start_position"], 
                              genes[gene_idx, "end_position"])) %>% 
    .$sample %>%
    unique()
  
  # Define region of interest around the gene of interest
  roi_margin = 4000000
  roi_gr = GRanges(
    seqnames = genes[gene_idx, "chromosome_name"],
    ranges = IRanges(
      start = genes[gene_idx, "start_position"] - roi_margin,
      end = genes[gene_idx, "end_position"] + roi_margin
    )
  )
  
  # Plot the genes within the region of interest
  genes_plotting_df =  getBM(attributes=c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "gene_biotype", "ensembl_gene_id"),
                             filters= "chromosome_name", 
                             values= gsub("chr", "", seqnames(roi_gr)), 
                             mart=ensembl)
  genes_plotting_df$my_name = ifelse(genes_plotting_df$hgnc_symbol != "", 
                                     genes_plotting_df$hgnc_symbol, 
                                     paste0(genes_plotting_df$ensembl_gene_id))
  genes_plotting_df = genes_plotting_df %>%
    filter(gene_biotype == "protein_coding") %>% 
    filter(hasOverlap_withChr(gsub("chr", "", as.character(seqnames(roi_gr))), start(roi_gr), end(roi_gr),
                              chromosome_name, start_position, end_position))

  genes_plotting_df$hgnc_symbol = 
    ifelse(genes_plotting_df$hgnc_symbol %in% c( 
      "FAM84B",
      "POU5F1B", 
      "MYC", 
      "TMEM75"), 
      genes_plotting_df$hgnc_symbol, 
      NA)
  
  genes.fig = genes_plotting_df %>%
    ggplot(aes(x=start_position, y=1)) +
    geom_rect(xmin=genes_plotting_df$start_position, xmax=genes_plotting_df$end_position, ymin=-Inf, ymax=0.1, color=NA, fill="black", alpha=0.5) + 
    ylim(0,2) +
    geom_text_repel(aes(x=start_position, y=0.1, label = hgnc_symbol), nudge_y=0.1, size=2, segment.size=0.1, min.segment.length = 0) +
    theme_kons2() + 
    theme(axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank()) +
    xlim(start(roi_gr),end(roi_gr))

  
  
  for (entity_idx in 1:length(entities)){
    
    this_seg = 
      seg %>% 
      filter(sample %in% ampl_samples) %>%
      filter(cancer_type == entities[entity_idx])
    
    # View(this_seg)
    # readline(prompt="Press [enter] to continue")
    
    if (nrow(this_seg) == 0) next
    
    n_ampl = length(unique(this_seg$sample))
    
    if (n_ampl < 5) next

    bins = tileGenome(seqlengths = seqlengths(BSgenome.Hsapiens.UCSC.hg19), 
                      tilewidth = 10000,
                      cut.last.tile.in.chrom = T)
    seqlevelsStyle(bins) = "NCBI"
    
    this_seg_grl = makeGRangesListFromDataFrame(this_seg, 
                             split.field = "sample",
                             start.field = "startpos",
                             end.field = "endpos",
                             keep.extra.columns = T)
    seqlevelsStyle(this_seg_grl) = "NCBI"
    seqlevels(bins, pruning.mode = "coarse") = seqlevels(this_seg_grl)
    
    seg_amp_cov = lapply(this_seg_grl, coverage)
    seg_amp_binned = lapply(seg_amp_cov, function (this_profile_cov) binnedAverage(bins, numvar = this_profile_cov, varname = "isAmp"))
    for (i in 1:length(seg_amp_binned)){
      seg_amp_binned[[i]]$Name = names(seg_amp_cov)[i]
    }
    seg_amp_binned_tb = do.call(rbind, lapply(seg_amp_binned, as_tibble))
    seg_amp_binned_tb$isAmp = seg_amp_binned_tb$isAmp > 0 # allow only 0 or 1 as bin values = all bins with amplicon overlap --> 1
    
    profile = seg_amp_binned_tb %>%
      group_by(seqnames, start, end) %>%
      summarise(n = sum(isAmp)) %>%
      ungroup()

    ### faster implementation but: adds incompletely covered bins ---> plots *slightly* incorrect
    # this_seg_cov = coverage(
    #   makeGRangesFromDataFrame(this_seg, 
    #                            start.field = "startpos",
    #                            end.field = "endpos",
    #                            keep.extra.columns = T
    #   ))
    # seqlevelsStyle(this_seg_cov) = "NCBI"
    # seqlevels(bins, pruning.mode = "coarse") = seqlevels(this_seg_cov)
    # 
    # profile = binnedAverage(bins, numvar = this_seg_cov, varname = "n") %>%
    #   as_tibble()
    
    f = profile %>%
      mutate(seqnames = as.character(seqnames)) %>% 
      filter(seqnames == genes[gene_idx, "chromosome_name"],
             start >= start(roi_gr), 
             end <= end(roi_gr)) %>%
      mutate(PercentAmplified = 100 * n / n_ampl) %>% 
      ggplot(aes(x=start, y=PercentAmplified)) + 
      geom_line() + 
      annotate(geom = "rect", 
               xmin=genes[gene_idx, "start_position"],
               xmax=genes[gene_idx, "end_position"],
               ymin=0, ymax=Inf) + 
      theme_kons2() + 
      xlim(start(roi_gr),end(roi_gr)) +
      ylab("Samples [%]") + 
      xlab(genes[gene_idx, "chromosome_name"])
    
    # print(f +
    #         ggtitle(paste0(genes[gene_idx, "hgnc_symbol"],
    #                        "\n",
    #                        "TCGA-", entities[entity_idx], " (n=", as.character(n_ampl), ")")))
    # readline(prompt="Press [enter] to continue")
    
    ggsave(
      paste0("/Volumes/Elements/MYCNAmplicon/Results/TCGA_ASCAT_Figs/",
             entities[entity_idx], "_", genes[gene_idx, "hgnc_symbol"], ".pdf"),
      egg::ggarrange(
        genes.fig +
          ggtitle(paste0(genes[gene_idx, "hgnc_symbol"],
                         "\n",
                         "TCGA-", entities[entity_idx], " (n=", as.character(n_ampl), ")")) + 
          theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x=element_blank()), 
        f,
        nrow=2,
        heights=c(0.5,2)
      ),
      height = 3, width = 3, useDingbats = F, onefile = F
    )
    
    this_seg %>%
      filter(chr == genes[gene_idx, "chromosome_name"]) %>%
      ggplot(aes(x=startpos, y=sample)) + 
      geom_errorbarh(aes(xmin=startpos, xmax=endpos), height=0) +
      annotate(geom = "rect", 
               xmin=genes[gene_idx, "start_position"],
               xmax=genes[gene_idx, "end_position"],
               ymin=0, ymax=Inf) +
      theme_kons2() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y=element_blank()) +
      ylab("") +
      xlab(paste0("chromosome ", genes[gene_idx, "chromosome_name"])) +
      ggtitle(paste0(genes[gene_idx, "hgnc_symbol"],
                     "\n",
                     "TCGA-", entities[entity_idx], " (n=", as.character(n_ampl), ")")) +
      facet_zoom(xlim = c(start(roi_gr), end(roi_gr))) +
      theme(strip.background = element_rect(fill="grey90", color=NA)) +
      ggsave(paste0("/Volumes/Elements/MYCNAmplicon/Results/TCGA_ASCAT_Figs/",
                    entities[entity_idx], "_", genes[gene_idx, "hgnc_symbol"], "_allamplicons.pdf"),
             height=2.5, width=3, useDingbats=F, onefile=F)
      
  }
}

# ------------------------------------------------------------------------------
# Chromosome 12 in GBM
# ------------------------------------------------------------------------------  

seg %>%
  filter(cancer_type == "GBM", chr == "12") %>% 
  ggplot(aes(x=startpos, y=sample)) + 
  geom_errorbarh(aes(xmin=startpos, xmax=endpos), height=0) +
  theme_kons2() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y=element_blank()) +
  ylab("") +
  ggtitle(paste0(genes[gene_idx, "hgnc_symbol"],
                 "\n",
                 entities[entity_idx]))

seg %>%
  filter(cancer_type == "GBM", chr == "12") %>%
  group_by(sample) %>%
  summarise(nFragments = n_distinct(startpos)) %>%
  ungroup() %>%
  group_by(nFragments) %>%
  summarise(nSamples = n_distinct(sample)) %>%
  ungroup() %>% 
  ggplot(aes(x=nFragments, y=nSamples)) + 
  geom_col()

single_cdk4_fragment_gbm = 
  seg %>%
  filter(cancer_type == "GBM") %>%
  mutate(CDK4amp = hasOverlap_withChr(chr, startpos, endpos, 12, 58141510, 58149796)) %>%
  group_by(sample) %>%
  summarise(nFragments = n_distinct(startpos),
            CDK4amp = any(CDK4amp)) %>%
  ungroup() %>%
  filter(nFragments == 1, CDK4amp == TRUE) %>%
  .$sample %>%
  unique()
length(single_cdk4_fragment_gbm)

# ------------------------------------------------------------------------------
# Chromosome 12 in SARC
# ------------------------------------------------------------------------------  

seg %>%
  filter(cancer_type == "SARC", chr == "12") %>% 
  ggplot(aes(x=startpos, y=sample)) + 
  geom_errorbarh(aes(xmin=startpos, xmax=endpos), height=0) +
  theme_kons2() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y=element_blank()) +
  ylab("") +
  annotate(geom = "rect", 
           xmin=genes[genes$hgnc_symbol == "CDK4", "start_position"],
           xmax=genes[genes$hgnc_symbol == "CDK4", "end_position"],
           ymin=0, ymax=Inf, color="red", alpha=0.5) +
  annotate(geom = "rect", 
           xmin=genes[genes$hgnc_symbol == "MDM2", "start_position"],
           xmax=genes[genes$hgnc_symbol == "MDM2", "end_position"],
           ymin=0, ymax=Inf, color="blue", alpha=0.5) +
  annotate(geom = "rect", 
           xmin=genes[genes$hgnc_symbol == "FRS2", "start_position"],
           xmax=genes[genes$hgnc_symbol == "FRS2", "end_position"],
           ymin=0, ymax=Inf) +
  annotate(geom = "rect", 
           xmin=genes[genes$hgnc_symbol == "HMGA2", "start_position"],
           xmax=genes[genes$hgnc_symbol == "HMGA2", "end_position"],
           ymin=0, ymax=Inf) +
  ggtitle("SARC chr12 amplifications")

# ------------------------------------------------------------------------------
# Investigate some examples without plateau
# ------------------------------------------------------------------------------
# 1       CCND1              11       69455855     69469242
# 2        CDK4              12       58141510     58149796
# 3        EGFR               7       55086714     55324313
# 4        MDM2              12       69201956     69239214
# 5         MYC               8      128747680    128753674
# 6        MYCN               2       16080686     16087129
# 7        TERC               3      169482308    169482848
# 8        TERT               5        1253262      1295184

gene_of_interest = "PDGFRA"
cancer_of_interest = "GBM"

single_amp_fragments = 
  seg %>%
  filter(cancer_type == cancer_of_interest) %>%
  mutate(targetamp = hasOverlap_withChr(chr, startpos, endpos, genes[genes$hgnc_symbol == gene_of_interest, "chromosome_name"], genes[genes$hgnc_symbol == gene_of_interest, "start_position"], genes[genes$hgnc_symbol == gene_of_interest, "end_position"])) %>%
  group_by(sample) %>%
  summarise(nFragments = n_distinct(startpos),
            targetamp = any(targetamp)) %>%
  ungroup() %>%
  filter(nFragments == 1, targetamp == TRUE) %>%
  .$sample %>%
  unique()
length(single_amp_fragments)

this_seg  = 
  seg %>%
  filter(sample %in% single_amp_fragments)

n_ampl = length(unique(this_seg$sample))

bins = tileGenome(seqlengths = seqlengths(BSgenome.Hsapiens.UCSC.hg19), 
                  tilewidth = 10000,
                  cut.last.tile.in.chrom = T)
seqlevelsStyle(bins) = "NCBI"

this_seg_cov = coverage(
  makeGRangesFromDataFrame(this_seg, 
                           start.field = "startpos",
                           end.field = "endpos",
                           keep.extra.columns = T
  ))
seqlevelsStyle(this_seg_cov) = "NCBI"
seqlevels(bins, pruning.mode = "coarse") = seqlevels(this_seg_cov)

profile = binnedAverage(bins, numvar = this_seg_cov, varname = "n") %>%
  as_tibble()

roi_margin = 1000000
profile %>%
  filter(seqnames == genes[genes$hgnc_symbol == gene_of_interest, "chromosome_name"]) %>%
  ggplot(aes(x=start, y=n)) + 
  geom_line() + 
  annotate(geom = "rect", 
           xmin=genes[genes$hgnc_symbol == gene_of_interest, "start_position"],
           xmax=genes[genes$hgnc_symbol == gene_of_interest, "end_position"],
           ymin=0, ymax=Inf) +
  xlim(genes[genes$hgnc_symbol == gene_of_interest, "start_position"]-roi_margin, genes[genes$hgnc_symbol == gene_of_interest, "end_position"]+roi_margin)+
  theme_kons2() + 
  ggtitle(paste0(gene_of_interest,
                 "\n",
                 cancer_of_interest, " Single Fragment"))
