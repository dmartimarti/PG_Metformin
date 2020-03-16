# Enrichment analysis for KEIO and ASKA library screens from Sim and Jen 
# analysis directory:  "/Users/dmarti14/Documents/MRC_Postdoc/Projects/Simren/Metformin_resistance"

library(tidyverse)
library(here)
library(openxlsx)
library(goseq)


# session options
options(width = 220)

# read and clean list of total gene in E coli K12 (from EcoCyc)
gen.db = read.table('databases/All_genes_of_E._coli_K-12_substr._MG1655_-_GO.txt', sep = '\t', header = TRUE)
gen.db = gen.db %>% select(-Accession.1, -Product, -UniProt)
total_genes = as.character(gen.db[,1])
total_genes = total_genes[-44]

# read list of genes
keio = read.table('KEIO/biolog/gene_names_KEIO.txt')
keio = as.character(t(keio))

# do it to test if each gene name is included in the gene database
cosa = as.integer(keio %in% total_genes)
names(cosa) = keio

# generate a vector with gene names and if they are 'DE' or not
# in our case it's not DE, it's Metformin resistant
gene.vector = as.integer(total_genes %in% keio)
names(gene.vector) = total_genes



# look at the supported organisms from the package
supportedOrganisms()




nullp(gene.vector , 'E. coli K12', id = 'geneSymbol')



goseq(  method = Hypergeometric)















