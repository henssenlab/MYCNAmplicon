library(GenomicAlignments)
library(plyranges)
library(dplyr)
library(ggplot2)
library(data.table)
source("/Volumes/Elements/MYCNAmplicon/Code/CustomThemes.R")

# ------------------------------------------------------------------------------
# Plot IMR-5/75 de novo assembly
# ------------------------------------------------------------------------------

paf = read.table("/Volumes/Elements/IMR-5-75.ONT-flye-assembly-meta-genomesize1g/assembly.maptohg19.paf") %>% as_tibble()
colnames(paf) = c(
  "Query", "QueryLength", "QueryStart", "QueryEnd", "RelativeStrand", 
  "ReferenceSequence", "ReferenceSequenceLength", "ReferenceSequenceStart", "ReferenceSequenceEnd",
  "NumberOfResidueMatches", "AlignmentBlockLength", "MappingQuality", 
  "NoIdea1", "NoIdea2", "NoIdea3", "NoIdea4", "NoIdea5"
)
paf = paf %>% filter(Query == "contig_499")
paf$ReferenceStart = ifelse(paf$RelativeStrand == "+", paf$ReferenceSequenceStart, paf$ReferenceSequenceEnd)
paf$ReferenceEnd = ifelse(paf$RelativeStrand == "+", paf$ReferenceSequenceEnd, paf$ReferenceSequenceStart)

paf %>%
  filter(ReferenceSequence == "chr2") %>% 
  ggplot() + 
  geom_segment(aes(x=QueryStart, xend=QueryEnd, y=ReferenceSequenceStart, yend=ReferenceSequenceEnd)) + 
  theme_kons2() +
  xlab("De novo assembly") + 
  ylab("chromosome 2 (hg19)") +
  ggtitle("IMR-5/75") +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/IMR575_NanoporeDeNovoAssembly_contig499_maptochr2.pdf", 
       height = 2, width=1.5)

paf %>%
  filter(ReferenceSequence == "chr2") %>% 
  ggplot() + 
  geom_segment(aes(y=QueryStart, yend=QueryEnd, x=ReferenceStart, xend=ReferenceEnd)) + 
  theme_kons2() +
  ylab("De novo assembly") + 
  xlab("chromosome 2 (hg19)") +
  ggtitle("Kelly") +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/IMR575_NanoporeDeNovoAssembly_contig499_maptochr2_Layout2.pdf", 
         height = 2, width=1.5)


# ------------------------------------------------------------------------------
# Plot Kelly de novo assembly
# ------------------------------------------------------------------------------

# samtools view -h -o /Volumes/Elements/Kelly.ONT-flye-assembly-meta-genomesize1g/assembly.maptohg19.sam /Volumes/Elements/Kelly.ONT-flye-assembly-meta-genomesize1g/assembly.maptohg19.bam
# ~/Desktop/nano-wgs/utils/k8-darwin ~/Desktop/nano-wgs/utils/sam2paf.js /Volumes/Elements/Kelly.ONT-flye-assembly-meta-genomesize1g/assembly.maptohg19.sam > /Volumes/Elements/Kelly.ONT-flye-assembly-meta-genomesize1g/assembly.maptohg19.paf

paf = read.table("/Volumes/Elements/Kelly.ONT-flye-assembly-meta-genomesize1g/assembly.maptohg19.paf") %>% as_tibble()
colnames(paf) = c(
  "Query", "QueryLength", "QueryStart", "QueryEnd", "RelativeStrand", 
  "ReferenceSequence", "ReferenceSequenceLength", "ReferenceSequenceStart", "ReferenceSequenceEnd",
  "NumberOfResidueMatches", "AlignmentBlockLength", "MappingQuality", 
  "NoIdea1", "NoIdea2", "NoIdea3", "NoIdea4", "NoIdea5"
)
paf = paf %>% filter(Query == "contig_31")
paf$ReferenceStart = ifelse(paf$RelativeStrand == "+", paf$ReferenceSequenceStart, paf$ReferenceSequenceEnd)
paf$ReferenceEnd = ifelse(paf$RelativeStrand == "+", paf$ReferenceSequenceEnd, paf$ReferenceSequenceStart)
paf %>%
  filter(ReferenceSequence == "chr2") %>% 
  ggplot() + 
  geom_segment(aes(x=QueryStart, xend=QueryEnd, y=ReferenceStart, yend=ReferenceEnd)) + 
  theme_kons2() +
  xlab("De novo assembly") + 
  ylab("chromosome 2 (hg19)") +
  ggtitle("Kelly") +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/Kelly_NanoporeDeNovoAssembly_contig_31_maptochr2.pdf", 
         height = 2, width=1.5)

paf %>%
  filter(ReferenceSequence == "chr2") %>% 
  ggplot() + 
  annotate(ymin=16080683, ymax=16087129, xmin=-Inf, xmax=Inf, geom="rect", color = "steelblue") +
  annotate(ymin=16375277, ymax=16385344, xmin=-Inf, xmax=Inf, geom="rect", color = "firebrick3") +
  geom_segment(aes(x=QueryStart, xend=QueryEnd, y=ReferenceStart, yend=ReferenceEnd)) + 
  theme_kons2() +
  xlab("De novo assembly") + 
  ylab("chromosome 2 (hg19)") +
  ggtitle("Kelly") +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/Kelly_NanoporeDeNovoAssembly_contig_31_maptochr2_Annot.pdf", 
         height = 2, width=1.5)
  
paf %>%
  ggplot() + 
  geom_segment(aes(x=QueryStart, xend=QueryEnd, y=ReferenceStart, yend=ReferenceEnd)) + 
  theme_kons2() +
  xlab("De novo assembly") + 
  facet_grid(ReferenceSequence ~ .) +
  ggtitle("Kelly") +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/Kelly_NanoporeDeNovoAssembly_contig_31_maptogenome.pdf", 
         height = 1.5, width=2)

paf %>%
  filter(ReferenceSequence == "chr2") %>% 
  ggplot() + 
  geom_segment(aes(y=QueryStart, yend=QueryEnd, x=ReferenceStart, xend=ReferenceEnd)) + 
  theme_kons2() +
  ylab("De novo assembly") + 
  xlab("chromosome 2 (hg19)") +
  ggtitle("Kelly") +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/Kelly_NanoporeDeNovoAssembly_contig_31_maptochr2_Layout2.pdf", 
         height = 2, width=1.5)



# ------------------------------------------------------------------------------
# Plot CHP-212 de novo assembly
# ------------------------------------------------------------------------------

paf = fread("/Volumes/Elements/CHP-212.ONT-flye-assembly-meta-genomesize1g/assembly.maptohg19.paf",
            fill=T) %>% as_tibble()
colnames(paf) = c(
  "Query", "QueryLength", "QueryStart", "QueryEnd", "RelativeStrand", 
  "ReferenceSequence", "ReferenceSequenceLength", "ReferenceSequenceStart", "ReferenceSequenceEnd",
  "NumberOfResidueMatches", "AlignmentBlockLength", "MappingQuality", 
  "NoIdea1", "NoIdea2", "NoIdea3", "NoIdea4", "NoIdea5", "NoIdea6"
)
paf = paf %>% filter(Query == "contig_1695")
paf$ReferenceStart = ifelse(paf$RelativeStrand == "+", paf$ReferenceSequenceStart, paf$ReferenceSequenceEnd)
paf$ReferenceEnd = ifelse(paf$RelativeStrand == "+", paf$ReferenceSequenceEnd, paf$ReferenceSequenceStart)
paf %>%
  filter(ReferenceSequence == "chr2") %>% 
  ggplot() + 
  geom_segment(aes(x=QueryStart, xend=QueryEnd, y=ReferenceStart, yend=ReferenceEnd)) + 
  theme_kons2() +
  xlab("De novo assembly") + 
  ylab("chromosome 2 (hg19)") +
  ggtitle("CHP-212") +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/CHP212_NanoporeDeNovoAssembly_contig1695_maptochr2.pdf", 
         height = 2, width=1.5)

paf %>%
  filter(ReferenceSequence == "chr2") %>% 
  ggplot() + 
  geom_segment(aes(y=QueryStart, yend=QueryEnd, x=ReferenceStart, xend=ReferenceEnd)) + 
  theme_kons2() +
  ylab("De novo assembly") + 
  xlab("chromosome 2 (hg19)") +
  ggtitle("CHP-212") +
  ggsave("/Volumes/Elements/MYCNAmplicon/Results/CHP212_NanoporeDeNovoAssembly_contig1695_maptochr2_Layout2.pdf", 
         height = 2, width=1.5)

 paf %>%
  ggplot() + 
  geom_segment(aes(x=QueryStart, xend=QueryEnd, y=ReferenceStart, yend=ReferenceEnd)) + 
  theme_kons2() +
  xlab("De novo assembly") + 
  facet_grid(ReferenceSequence ~ .) +
  ggtitle("CHP-212")

 # ------------------------------------------------------------------------------
 # Plot NGP de novo assembly
 # ------------------------------------------------------------------------------
 
 paf = fread("/Volumes/Elements/NGP.ONT-flye-assembly-meta-genomesize1g/assembly.maptohg19.paf",
             fill=T) %>% as_tibble()
 colnames(paf) = c(
   "Query", "QueryLength", "QueryStart", "QueryEnd", "RelativeStrand", 
   "ReferenceSequence", "ReferenceSequenceLength", "ReferenceSequenceStart", "ReferenceSequenceEnd",
   "NumberOfResidueMatches", "AlignmentBlockLength", "MappingQuality", 
   "NoIdea1", "NoIdea2", "NoIdea3", "NoIdea4", "NoIdea5", "NoIdea6"
 )
 paf = paf %>% filter(Query == "contig_45")
 paf$ReferenceStart = ifelse(paf$RelativeStrand == "+", paf$ReferenceSequenceStart, paf$ReferenceSequenceEnd)
 paf$ReferenceEnd = ifelse(paf$RelativeStrand == "+", paf$ReferenceSequenceEnd, paf$ReferenceSequenceStart)
 
 paf %>%
   filter(ReferenceSequence == "chr2") %>% 
   ggplot() + 
   annotate(ymin=16080683, ymax=16087129, xmin=-Inf, xmax=Inf, geom="rect", color = "steelblue") +
   annotate(ymin=16375277, ymax=16385344, xmin=-Inf, xmax=Inf, geom="rect", color = "firebrick3") +
   geom_segment(aes(x=QueryStart, xend=QueryEnd, y=ReferenceStart, yend=ReferenceEnd)) + 
   theme_kons2() +
   xlab("De novo assembly") + 
   ylab("chromosome 2 (hg19)") +
   ggtitle("NGP") +
   ylim(15000000,18000000) +
   ggsave("/Volumes/Elements/MYCNAmplicon/Results/NGP_NanoporeDeNovoAssembly_contig45_maptochr2_annotateMYCNAndE4.pdf", 
          height = 2, width=1.5)
 
 paf %>%
   filter(ReferenceSequence == "chr2") %>% 
   ggplot() + 
   geom_segment(aes(x=QueryStart, xend=QueryEnd, y=ReferenceStart, yend=ReferenceEnd)) + 
   theme_kons2() +
   xlab("De novo assembly") + 
   ylab("chromosome 2 (hg19)") +
   ggtitle("NGP") +
   ylim(15000000,18000000) +
   ggsave("/Volumes/Elements/MYCNAmplicon/Results/NGP_NanoporeDeNovoAssembly_contig45_maptochr2.pdf", 
          height = 2, width=1.5)
 
 paf %>%
   filter(ReferenceSequence == "chr2") %>% 
   ggplot() + 
   geom_segment(aes(y=QueryStart, yend=QueryEnd, x=ReferenceStart, xend=ReferenceEnd)) + 
   theme_kons2() +
   ylab("De novo assembly") + 
   xlab("chromosome 2 (hg19)") +
   ggtitle("NGP") +
   ggsave("/Volumes/Elements/MYCNAmplicon/Results/NGP_NanoporeDeNovoAssembly_contig45_maptochr2_Layout2.pdf", 
          height = 2, width=1.5)
 
 paf %>%
   ggplot() + 
   geom_segment(aes(x=QueryStart, xend=QueryEnd, y=ReferenceStart, yend=ReferenceEnd)) + 
   theme_kons2() +
   xlab("De novo assembly") + 
   facet_grid(ReferenceSequence ~ .) +
   ggtitle("NGP") +
   ggsave("/Volumes/Elements/MYCNAmplicon/Results/NGP_NanoporeDeNovoAssembly_contig45_allchr.pdf", 
          height = 2, width=1.5)
 
 
 
