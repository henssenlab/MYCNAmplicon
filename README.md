Data analysis for Helmsauer, Valieva, Ali et al. *Enhancer hijacking determines extrachromosomal MYCN amplicon architecture in neuroblastoma* ([bioRxiv](https://www.henssenlab.com))

- `IdentifyCandidateEnhancers.R` integrates RNA-seq and ChIP-seq to identify *MYCN*-driving enhancers in neuroblastoma
- `DepuydtDataMYCNAmplicon.R` generates aggregate copy number profiles around *MYCN*
- `AmpliconEnhancerCharacterization.R` plots epigenetic data on *MYCN* amplicon.
- `ChIPSeqTracksOverviewFigures.R` plots H3K27ac ChIP-seq data around *MYCN* for different neuroblastoma cell lines
- `RadaIglesias.R` plots H3K27ac ChIP-seq around *MYCN* for developmental cell types
- `MotifScanning.R` scans enhancers for transcription factor binding motifs
- `ReconstructAmplicon_{CHP212,IMR575}.R` reconstructs the *MYCN* amplicon from structural rearrangement calls for neuroblastoma cell lines and plots the reconstruction together with epigenomic data
- `FragileSites.R` overlaps aggregate copy number profiles and common fragile sites
- `ClassIIWGSRearrangementPlotting.R` plots rearrangement data for three patients that lack the e4 enhancer
- `EnhancerVsCRCChIPseq.R` describes overlap between core regulatory circuit factor ChIP-seq and *MYCN*-driving enhancers
- `prepareRNAseq.R` parses and normalizes RNA-seq data
- `PlotAssemblyMapping.R` plots mapping of *de novo* assembled contigs to hg19

Please do not hesitate to contact us at henssenlab@gmail.com.
