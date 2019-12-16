library(dplyr)
library(parallel)
cores = detectCores()
source("/Volumes/Elements/MYCNAmplicon/Code/hasOverlap.R")

# ------------------------------------------------------------------------------
# Parse VCF files from Svaba Tumor-Normal
# ------------------------------------------------------------------------------

parse_svaba_allsv = function(indels_fname, svfromindels_fname, sv_fname, sample=NA, cohort=NA, determine_type=FALSE){

  # indels_fname = NULL
  # svfromindels_fname = NULL
  # sv_fname = "/Volumes/Elements/MYCNAmplicon/Data/Konstantin_unfiltered_svaba_results_21.10.19/NB2013/somatic_run.svaba.unfiltered.somatic.sv.vcf"
  # sample = "NB2013"
  # cohort = NA
  # determine_type = F
  
  # indels_fname = NULL
  # svfromindels_fname = NULL
  # sv_fname = "/Volumes/Elements/MYCNAmplicon/Data/Konstantin_unfiltered_svaba_results_21.10.19/CB2003/somatic_run.svaba.unfiltered.somatic.sv.vcf"
  # sample = "CB2003"
  # cohort = NA
  # determine_type = F

  indels_vcf = tryCatch(read.table(indels_fname,
                                   comment.char = "#",
                                   header = F,
                                   sep = "\t") %>% as_tibble(),
                        error = function(e) NULL)
  
  if (!is.null(indels_vcf)){
    colnames(indels_vcf) = c("ChrA", "PosA", "ID","REF","ALT","QUAL","FILTER", "INFO", "FORMAT", "NORMAL", "TUMOR")
    indels_vcf$Sample = sample
    indels_vcf$Cohort = cohort
    indels_vcf$Class = "indel"
    indels_vcf$ChrA = as.character(indels_vcf$ChrA)
    indels_vcf$PosA = as.numeric(indels_vcf$PosA)
    indels_vcf$ChrB = indels_vcf$ChrA
    indels_vcf$PosB = indels_vcf$PosA + nchar(as.character(indels_vcf$REF)) - 1 # checked, those are exactly the ucsc
    indels_vcf$Filter = indels_vcf$FILTER
    
    indels_vcf = indels_vcf %>% 
      dplyr::select(Sample, Cohort, Class, ChrA, PosA, ChrB, PosB, Filter)
  }
  
  ############################
  
  svfromindels_vcf = tryCatch(read.table(svfromindels_fname,
                                         comment.char = "#",
                                         header = F,
                                         sep = "\t") %>% as_tibble(),
                              error = function(e) NULL)
  
  if (!is.null(svfromindels_vcf)){
    colnames(svfromindels_vcf) = c("ChrA",
                                   "PosA",
                                   "ID",
                                   "REF",
                                   "ALT",
                                   "QUAL",
                                   "FILTER",
                                   "INFO",
                                   "FORMAT",
                                   "NORMAL", 
                                   "TUMOR")
    svfromindels_vcf$Sample = sample
    svfromindels_vcf$Cohort = cohort
    svfromindels_vcf$Class = "indelfromsv"
    svfromindels_vcf$ChrA = as.character(svfromindels_vcf$ChrA)
    svfromindels_vcf$PosA = as.numeric(svfromindels_vcf$PosA)
    svfromindels_vcf$ChrB = svfromindels_vcf$ChrA
    svfromindels_vcf$PosB = svfromindels_vcf$PosA + nchar(as.character(svfromindels_vcf$REF)) - 1
    svfromindels_vcf = svfromindels_vcf %>% 
      mutate(Filter=FILTER) %>%
      dplyr::select(Sample, Cohort, Class, ChrA, PosA, ChrB, PosB, Filter)
  }
  
  ############################
  
  sv_vcf = tryCatch(read.table(sv_fname,
                               comment.char = "#",
                               header = F,
                               sep = "\t") %>% as_tibble(),
                    error = function(e) NULL)

  if (!is.null(sv_vcf)){
    
    colnames(sv_vcf) = c("ChrA","PosA","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","NORMAL", "TUMOR")
    
    sv_vcf$BNDPairID = sv_vcf$ID %>% gsub(":.*$", "", .)
    sv_vcf$ChrA = as.character(sv_vcf$ChrA)
    sv_vcf$PosA = as.numeric(sv_vcf$PosA)
    sv_vcf$Sample = sample
    
    if (determine_type){
      sv_vcf$Type = sapply(sv_vcf$ID, function(ID) svaba_get_sv_type(ID, sv_vcf))
    } else {
      sv_vcf$Type = "SV"
    }
    sv_vcf$BracketCoord = sv_vcf$ALT
    sv_vcf$DirectionA = ifelse(grepl("^[AGCTN]", sv_vcf$ALT), "Tail", "Head") # Head = Breakpoint Coordinate is Max Coordinate = Coordinate decrease with distance from Breakpoint
    sv_vcf$DirectionB = ifelse(grepl("\\[", sv_vcf$ALT), "Head", "Tail")
    sv_vcf$JunctionType = paste0(sv_vcf$DirectionA, ">", sv_vcf$DirectionB)

    sv_vcf$ALT = ifelse(
      grepl("\\[", sv_vcf$ALT),
      sv_vcf$ALT %>% gsub("\\[[[:alnum:]]*$", "", .) %>% gsub("^[[:alnum:]]*\\[", "", .),
      sv_vcf$ALT %>% gsub("\\][[:alnum:]]*$", "", .) %>% gsub("^[[:alnum:]]*\\]", "", .)
    )
    
    sv_vcf$ChrB = gsub(":.*$", "", sv_vcf$ALT) %>% as.character() # some go wrong, because
    sv_vcf$PosB = gsub("^.*:", "", sv_vcf$ALT) %>% as.numeric() # some go wrong,
    sv_vcf$Evidence = sv_vcf$INFO %>% gsub("^.*EVDNC=", "", .) %>% gsub(";.*", "", .) %>% as.character()
    sv_vcf$Homology = ifelse(
      grepl("HOMSEQ", sv_vcf$INFO),
      sv_vcf$INFO %>% gsub("^.*;HOMSEQ=", "", .) %>% gsub(";.*", "", .) %>% as.character(),
      ""
    )
    sv_vcf$Homology = ifelse(grepl(";IMPRECISE", sv_vcf$INFO),
                             NA,
                             sv_vcf$Homology)
    sv_vcf$HomologyLength = nchar(sv_vcf$Homology)

    sv_vcf$Insertion = ifelse(
      grepl("INSERTION", sv_vcf$INFO),
      sv_vcf$INFO %>% gsub("^.*;INSERTION=", "", .) %>% gsub(";.*", "", .) %>% as.character(),
      ""
    )
    sv_vcf$Insertion = ifelse(grepl(";IMPRECISE", sv_vcf$INFO),
                             NA,
                             sv_vcf$Insertion)
    sv_vcf$InsertionLength = nchar(sv_vcf$Insertion)

    sv_vcf$Class = "sv"
    sv_vcf$Cohort = cohort
    #sv_vcf$Type = sv_vcf$INFO %>% gsub("^.*;SVTYPE=", "", .) %>% gsub(";.*", "", .) %>% as.character()
    sv_vcf$MappingQuality = sv_vcf$INFO %>% gsub("^.*;MAPQ=", "", .) %>% gsub(";.*", "", .) %>% as.numeric()
    sv_vcf$Quality = sv_vcf$QUAL %>% as.numeric()
    
    if (length(unique(sv_vcf$FORMAT)) > 1){
      stop("FORMAT row is not identical for different SVs.")
    } else {
      sv_vcf_format = as_tibble(do.call(rbind, strsplit(as.character(sv_vcf$FORMAT), "\\:", perl=T)))
      sv_vcf_normal = as_tibble(do.call(rbind, strsplit(as.character(sv_vcf$NORMAL), "\\:", perl=T)))
      colnames(sv_vcf_normal) = paste0("Normal_",sv_vcf_format[1,])
      sv_vcf_tumor = as_tibble(do.call(rbind, strsplit(as.character(sv_vcf$TUMOR), "\\:", perl=T)))
      colnames(sv_vcf_tumor) = paste0("Tumor_",sv_vcf_format[1,])
    }
    sv_vcf = cbind(sv_vcf, sv_vcf_normal, sv_vcf_tumor)
    
    sv_vcf = 
      sv_vcf %>%
      filter(grepl(":1$", ID)) # assuming that *:1 for intrachromosomal rearrangements always is PosA < PosB
    
    sv_vcf = sv_vcf %>% mutate(Filter=FILTER, Tumor_LocalCoverage=as.numeric(Tumor_DP), TumorSplitReads=as.numeric(Tumor_SR), TumorDiscordantReads=as.numeric(Tumor_DR), TumorReadSupport=as.numeric(Tumor_AD),
                               NormalSplitReads=as.numeric(Normal_SR), NormalDiscordantReads=as.numeric(Normal_DR), NormalReadSupport=as.numeric(Normal_AD)) %>%
      dplyr::select(Sample, Cohort, Class, BNDPairID, ID, ChrA, PosA, ChrB, PosB, Filter, Type, DirectionA, DirectionB, JunctionType, Evidence, Tumor_LocalCoverage, TumorReadSupport, TumorDiscordantReads, TumorSplitReads, NormalReadSupport, NormalDiscordantReads, NormalSplitReads, Homology,HomologyLength,Insertion,InsertionLength)
  }
  
  vcf = bind_rows(sv_vcf, indels_vcf, svfromindels_vcf) %>% as_tibble()
  
  vcf$ChrA  = vcf$ChrA %>% gsub("chr", "", .) %>% paste0("chr", .)
  vcf$ChrB  = vcf$ChrB %>% gsub("chr", "", .) %>% paste0("chr", .)

  return(vcf)
}

svaba_get_sv_type <- function(svaba_vcf_table_ID, svaba_vcf_table){
  # adapted from: https://github.com/walaj/svaba/issues/4
  
  #svaba_vcf_table_ID = sv_vcf[[1, "ID"]]
  #svaba_vcf_table = sv_vcf
  
  root <- gsub(":[12]", "", svaba_vcf_table_ID)
  mate1 <- paste0(root, ":1")
  mate2 <- paste0(root, ":2")
  
  alt1 = tryCatch(svaba_vcf_table %>% filter(ID == mate1) %>% .$ALT, error = function(e) NA) %>% as.character()
  alt2 = tryCatch(svaba_vcf_table %>% filter(ID == mate2) %>% .$ALT, error = function(e) NA) %>% as.character()
  chr1 = tryCatch(svaba_vcf_table %>% filter(ID == mate1) %>% .$ChrA, error = function(e) NA) %>% as.character()
  chr2 = tryCatch(svaba_vcf_table %>% filter(ID == mate2) %>% .$ChrA, error = function(e) NA) %>% as.character()
  
  if (identical(alt1, character(0)) | identical(alt2, character(0)) | identical(chr1, character(0)) | identical(chr2, character(0))) return(NA)
  if (is.na(alt1) | is.na(alt2) | is.na(chr1) | is.na(chr2)) return(NA)
  
  if ((grepl("\\[", alt1) & grepl("\\[", alt2)) | (grepl("\\]", alt1) & grepl("\\]", alt2))){
    sv_type <- "INV"
    
  } else if (grepl("[ACTGN]\\[", alt1) & grepl("^\\]", alt2)){
    sv_type <- "DEL"
    
  } else if (grepl("^\\]", alt1) & grepl("[ACTGN]\\[", alt2)){
    sv_type <- "DUP/INS"
    
  } else{
    sv_type <- "UNKNOWN"
  }
  
  # own category BND for interchromosomal rearrangements
  if (chr1 != chr2) sv_type = "BND"
  
  return(sv_type)
}

# ------------------------------------------------------------------------------
# Read all (unfiltered) svaba calls
# ------------------------------------------------------------------------------

datapath = "/Volumes/Elements/MYCNAmplicon/Data/Konstantin_unfiltered_svaba_results_21.10.19/"
samples = dir(datapath)
samples = samples[!grepl("excluded", samples)]
#samples = samples[!grepl("relapse", samples)] # CB2034_relapse is *huge* (>20fold larger than the usual vcf file)
#samples = samples[samples != "NBL50"]
#samples = samples[samples != "NBL55"]

svaba = list()
for (i in 1:length(samples)){
  print(samples[i])
  svaba[[i]] = parse_svaba_allsv(NULL, NULL, paste0(datapath, samples[i], "/somatic_run.svaba.unfiltered.somatic.sv.vcf"), samples[i])
}
svaba = do.call(rbind, svaba) 
svaba$Sample = gsub("CB", "NB", svaba$Sample)
samples = gsub("CB", "NB", samples)

# ------------------------------------------------------------------------------
# Add Blacklist Information
# ------------------------------------------------------------------------------

blacklist = read.table("/Volumes/Elements/MYCNAmplicon/Data/hg19-blacklist.v2.bed",
                       header = F, sep="\t")
colnames(blacklist) = c("Chr", "Start", "End", "Class")
isblacklisted = function(chr, pos) any(hasOverlap_withChr(chr, pos, pos, blacklist$Chr, blacklist$Start, blacklist$End))
svaba$isBlacklistedA = mcmapply(isblacklisted, svaba$ChrA, svaba$PosA, mc.cores=cores)
svaba$isBlacklistedB = mcmapply(isblacklisted, svaba$ChrB, svaba$PosB, mc.cores=cores)

# ------------------------------------------------------------------------------
# Save Rdata
# ------------------------------------------------------------------------------

rm(list = setdiff(ls(), c("samples", "svaba")))
save.image("/Volumes/Elements/MYCNAmplicon/Data/UnfilteredSvabaSVs.Rdata")

