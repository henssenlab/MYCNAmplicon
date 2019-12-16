Code accompanying the Helmsauer, Valieva, Ali et al. *Enhancer hijacking determines extrachromosomal MYCN amplicon architecture in neuroblastoma* ([bioRxiv](https://www.henssenlab.com))

- `IdentifyCandidateEnhancers.R` integrates RNA-seq and ChIP-seq to identify *MYCN*-driving enhancers in neuroblastoma
- `DepuydtDataMYCNAmplicon.R` generates aggregate copy number profiles around *MYCN*
- `AmpliconEnhancerCharacterization.R` plots epigenetic data on *MYCN* amplicon.
- `ChIPSeqTracksOverviewFigures.R` plots H3K27ac ChIP-seq data around *MYCN* for different neuroblastoma cell lines
- `RadaIglesias.R` plots H3K27ac ChIP-seq around *MYCN* for developmental cell types
- `MotifScanning.R` scans enhancers for transcription factor binding motifs
- `ReconstructAmplicon_{CHP212,IMR575}.R` reconstructs the *MYCN* amplicon from structural rearrangement calls for neuroblastoma cell lines and plots the reconstruction together with epigenomic data
- `FragileSites.R` overlaps aggregate copy number profiles and common fragile sites

Please do not hesitate to contact us at henssenlab@gmail.com.
