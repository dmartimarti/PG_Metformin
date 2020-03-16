## network printing from KEIO and ASKA gene enrichment

# setwd('/Users/dmarti14/Documents/MRC_Postdoc/Projects/Simren/Metformin_resistance/gene_network')

# load libraries
library(tidyverse)
library(readxl)
library(igraph)

options(width = 240)

# read files



keio.genes = read.table('gene_names_KEIO.txt', sep = '\t', header = FALSE) # list of KEIO genes
aska.genes = read.table('gene_names_ASKA.txt', sep = '\t', header = FALSE) # list of ASKA genes
keio.net = read.table('keio_net.tsv', sep = '\t', header = TRUE) # keio connections
aska.net = read.table('aska_net.tsv', sep = '\t', header = TRUE) # aska connections
all.net  =  read.table('all_net.tsv', sep = '\t', header = TRUE) # all genes connections
coord = read.delim('all_network_coordinates.txt', sep = '\t', header = TRUE, quote = "") # coordinates of all genes

keio.genes = keio.genes %>%
	mutate(set = 'KEIO', 
		   setn = 1,
		   V1 = as.character(V1)) 

aska.genes = aska.genes %>%
	mutate(set = 'ASKA',
		   setn = 2,
		   V1 = as.character(V1))


all.genes = unique(rbind(keio.genes, aska.genes))
names(all.genes) = c('genes', 'set', 'setn')
# # list of genes
# keio.genes = as.character(keio.genes[,1])
# aska.genes = as.character(aska.genes[,1])

# all.genes = c(keio.genes, aska.genes)


# it seems no all genes in all.net are included in all.genes
all.net[,1] %in% all.genes[,1]

all.net[,1][!(all.net[,1] %in% all.genes[,1])]
all.net[,2][!(all.net[,2] %in% all.genes[,1])]

all.net = all.net %>% 
	mutate(node1 = as.character(node1),
		   node2 = as.character(node2))

# fixing gene names (because of the synonyms)
all.net[,1][all.net[,1] == 'hofN'] = 'yrfC'
all.net[,1][all.net[,1] == 'mnmA'] = 'trmU'
all.net[,1][all.net[,1] == 'rarA'] = 'ycaJ'
all.net[,2][all.net[,2] == 'hofN'] = 'yrfC'
all.net[,2][all.net[,2] == 'mnmA'] = 'trmU'
all.net[,2][all.net[,2] == 'rarA'] = 'ycaJ'

# build the igraph object
net = graph_from_data_frame(d = all.net, vertices = all.genes, directed = F) 


colrs <- c("gray50", "tomato")

V(net)$color = colrs[V(net)$setn]


# plot graph
plot(net, edge.arrow.size = .9,
		  vertex.size = 5,
		  # vertex.label.dist = 2,
		  vertex.label.cex = 0.7)




















