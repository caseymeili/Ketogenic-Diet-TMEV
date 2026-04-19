# Ketogenic Diet TMEV

### Amplicon Processing Pipeline and Alpha Diversity
Mothur-based pipeline (v1.48.3) for processing paired-end 16S rRNA sequencing data (V3–V4). Includes read assembly, quality filtering, alignment to SILVA, chimera removal (VSEARCH), taxonomic classification (RDP v19), OTU clustering, and generation of shared and taxonomy files for downstream analysis. Also performs subsampling, rarefaction, and alpha diversity calculations.

**Required File for Sequence Processings**
- SILVA database (release 138.2) obtained from: https://mothur.org/wiki/silva_reference_files/
- PDS reference files (version 19) obtained from: https://mothur.org/wiki/rdp_reference_files/


### Beta Diversity Analysis
OTU abundance tables and taxonomy files generated from mothur were compiled into a Phyloseq object (see make_phyloseq.R for construction of the OTU table), and sample metadata were incorporated. This script performs beta diversity analyses, including Bray-Curtis distance calculation, principal coordinates analysis (PCoA) ordination, PERMANOVA, and assessment of homogeneity of dispersion (betadisper).


### Differential Abundance Analysis


### Questions 
Questions can be directed to Casey Meili (casey.meili@pharm.utah.edu)


------------------------------------------------------------------------------------------------


**References**

Schloss, P.D., et al., Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Appl Environ Microbiol, 2009. 75(23):7537-41

Kozich JJ, Westcott SL, Baxter NT, Highlander SK, Schloss PD. (2013): Development of a dual-index sequencing strategy and curation pipeline for analyzing amplicon sequence data on the MiSeq Illumina sequencing platform. Applied and Environmental Microbiology. 79(17):5112-20.


