# PG_Metformin
Git repository for data analysis of E. coli pangenomic library for Metformin resistance

Scripts are divided into two folders, with the main scripts listed here:

**phylo**: scripts created to work with the _E. coli_ genomes and annotations.  
- **files_fetch.py**: script to download all genomes and annotations from the original database
- **anvio_pan.sh**: basic instructions to create a small pan-genome with Anvio
- **genomic_distances.R**: to calculate and plot genomic distances from the Phylonium python package
- **tree_representation.R**: master script to represent the phylogenetic tree of the core genome from the pan-genomic library

**R_data_analysis**: scripts meant to analyse data non-related to the phylogenetic tree:
- **resistance_growth.R**: script to analyse PG 
- **chem_space_PCA.R**: chemical space analysis with Biolog compounds, also present [here](https://dmartimarti.github.io/Chem_space/) as a github-page
- **synergy_enrichment.R**: drug-drug synergy analysis from SynergyFinder results 
- **metf_40_60_biolog_clean.R** and **metf_40_60_biolog.R**: different versions of the analysis from the metformin-drug bacterial growth from Biospa and Biolog

Everything else still needs to be tidied, classified, and placed in other folders.
