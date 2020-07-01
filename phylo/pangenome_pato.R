### PATO workflow for pangenomes

library(pato)
library(tidyverse)

# load vignettes 
browseVignettes("pato")

## test 

setwd('/home/dani/Documents/MRC_postdoc/Pangenomic/phylo/original_data/pangn_ffn')

files <- read.table("my_list.txt",header = FALSE)

# create main files from mash and mmseqs2
ecoli_mash_all <- mash(files, n_cores = 12, type = 'nucl')
ecoli_mm<- mmseqs(files, coverage = 0.8, identity = 0.8, evalue = 1e-6, n_cores = 12)

ecoli_accnet_all <- accnet(ecoli_mm, threshold = 0.8, singles = FALSE)


# detect outliers (not the same species)
outl <- outliers(ecoli_mash_all,threshold = 0.06)

# if you have outliers, you can remove them with
ecoli_mash_all <-remove_outliers(ecoli_mash_all, outl)
ecoli_accnet_all <-remove_outliers(ecoli_accnet_all, outl)

# remove files that were outliers from your original list
files <-  anti_join(files, outl, by=c("V1"="Source"))


# For select 800 non redundant samples
nr_list <- non_redundant(ecoli_mash_all, distance = 0.0001)

# to create the objects only with the representatives of each cluster:
efaecium_accnet <- extract_non_redundant(efaecium_accnet_all, nr_list)
efaecium_mash <- extract_non_redundant(ecoli_mash_all, nr_list)


# Core genome analysis
cp <- core_plots(ecoli_mm,reps = 12, threshold = 0.95, steps = 10)

# calculate core genome by hand
ecoli_mm<- mmseqs(files, coverage = 0.8, identity = 0.8, evalue = 1e-6, n_cores = 12)

df = ecoli_mm$table %>% as_tibble()

nGenomes  =  df  %>% distinct(Genome_genome) %>% count()

table <-  df %>%
  distinct() %>%
  group_by(Prot_prot) %>%
  mutate(Nprot = n_distinct(Genome_prot), 
         Ngenomes = n_distinct(Genome_genome)) %>% 
  ungroup() %>%
  filter(Ngenomes > nGenomes$n*0.99)

table



core <- mmseqs(files, n_cores = 12, coverage = 0.8, identity = 0.8) 
core = core_genome(ecoli_mm, type = 'nucl', n_cores = 12)

export_core_to_fasta(core,"core.aln")

# we can read the core tree and plot it 
# iqtree -s core.aln -fast -nt 12 -m GTR+I+G

core_tree = ape::read.tree("core.aln.treefile")
core_tree %>% phytools::midpoint.root() %>% ggtree::ggtree()


# accessory genome
efaecium_accnet_all <- accnet(efaecium_mm,threshold = 0.8, singles = FALSE)

ef_nr<- non_redundant(efaecium_mash,distance = 0.0001)
efaecium_accnet_nr <- extract_non_redundant(efaecium_accnet, ef_nr)
ef_800_cl <- clustering(efaecium_accnet_nr, method = "mclust", d_reduction = TRUE)

export_to_gephi(efaecium_accnet_nr, "accnet800", cluster = ef_800_cl)

# enrichment analysis
accnet_enr_result <- accnet_enrichment_analysis(efaecium_accnet_nr, cluster = ef_800_cl)

accnet_enr_result

# Now, we can export a new network with the adjusted p-values as edge-weigth.
accnet_with_padj(accnet_enr_result) %>% export_to_gephi("accnet800.padj", cluster = ef_800_cl)


# knnn method of visualisation

# K-NNN with 10 neighbours and with repetitions
knnn_mash_10_w_r <- knnn(efaecium_mash,n=10, repeats = TRUE)
# K-NNN with 25 neighbours and with repetitions
knnn_mash_25_w_r <- knnn(efaecium_mash,n=25, repeats = TRUE)
# K-NNN with 50 neighbours and with repetitions
knnn_mash_50_w_r <- knnn(efaecium_mash,n=50, repeats = TRUE)

export_to_gephi(knnn_mash_50_w_r,file = "knnn_50_w_r.tsv")


ef_cl_mclust_umap <- clustering(efaecium_mash, method = "mclust",d_reduction = TRUE)
ef_cl_knnn <-clustering(knnn_mash_50_w_r, method = "louvain")
umap_mash <- umap_plot(efaecium_mash)
umap_mash <- umap_plot(efaecium_mash, cluster = ef_cl_mclust_umap)

cl_louvain = clustering(knnn_mash_25_w_r, method = "louvain")
plot_knnn_network(knnn_mash_25_w_r, cluster = cl_louvain, edge.alpha = 0.1)
