
# libraries ####

library(tidyverse)
library(tidytree)
library(treeio)
library(seqinr)
library(ggtree)
library(openxlsx)
library(treeio)
library(tidytree)
library(here)
library(readxl)
library(colorspace)
library(phytools)


# load data ---------------------------------------------------------------



# load straindb
strain_db = read_excel("D:/MRC_Postdoc/Pangenomic/phylo/original_data/strain_db.xlsx")

strain_db = strain_db %>% unite(names, ID, Assembly, remove = FALSE) 

strain_db2 = strain_db %>% 
  filter(Broadphenotype != 'Evolutionexperiment')

# load tree
tree = read.tree('core_gene_alignment.aln.treefile')

# summary statistics from worm imaging
worm.img = read_excel("D:/MRC_Postdoc/Pangenomic/Worm_imaging/summary/worm_imaging_stats.xlsx", 
                                 sheet = "Stats_per_strain")

# load AUC sum from PG library screen
resistance = read_excel("200601_2v3_top50_sensitivePG.xlsx", 
                                            sheet = "111")


# general tree ------------------------------------------------------------

# phylogroup info
sp_names = c(tree$tip.label)
df = data.frame(sp_names)
colnames(df) = 'label'

# check
strain_db %>% 
  filter(names %in% sp_names)

# tree %>% as_tibble() %>% full_join(df) %>% view()


# add phylogroup info and fix NT12139 strain
df = df %>% left_join(strain_db %>% 
                   select(names, phylogroup) %>%
                   rename(label = names)) %>%
  mutate(phylogroup = case_when(label == 'NT12139_282' ~ 'C',
                                label != 'NT12139_282'~ phylogroup))

# root tree and specify groups 
tree1 = tree %>% root('NT12461_14') %>% full_join(df) 

# plot general version of the tree
ggtree(tree1, layout = 'fan', size = 0.4, color = 'grey') +
  geom_tiplab(size = 1.5, aes(color = phylogroup, angle = angle, show.legend = FALSE)) +
  scale_colour_manual(values = c(
    "#E41A1C", # A
    "grey90",
    "#377EB8", # B1
    "#4DAF4A", # B2
    "#984EA3", # C
    "grey90",
    "grey90",
    "grey90",
    "grey90",
    "grey90",
    "#FF7F00", # D
    "#FFFF33", # E
    "#FFFF33", # E or clade I
    "#A65628", # F
    "grey90", 
    "#F781BF",  # G
    'grey'
  )) + guides(colour = guide_legend(override.aes = list(size = 6, pch = 18)))

ggsave(file = here('Rplots', 'rtree_phylogroups.pdf'), width = 100, height = 100, units = 'mm', scale = 2, device = 'pdf')



## remove tips from clades

to_drop = strain_db %>% 
  filter(phylogroup %in% c('cladeI', 'cladeII', 'cladeIII', 'cladeIV', 'cladeV', 'E or cladeI', 'albertii')) %>% 
  select(names) %>% t %>% as.character

nhx_reduced <- treeio::drop.tip(tree1, to_drop)

ggtree(nhx_reduced, layout = 'circular', size = 0.4, color = 'grey') +
  geom_tiplab(size = 1.5, aes(color = phylogroup, angle = angle, show.legend = FALSE)) +
  scale_colour_manual(values = c(
    "#E41A1C", # A
    "#377EB8", # B1
    "#4DAF4A", # B2
    "#984EA3", # C
    "#FF7F00", # D
    "#FFFF33", # E
    "#A65628", # F
    "grey90", 
    "#F781BF",  # G
    'grey'
  )) + guides(colour = guide_legend(override.aes = list(size = 6, pch = 18))) 

ggsave(file = here('Rplots', 'rtree_phylogroups_noClades.pdf'), width = 100, height = 100, units = 'mm', scale = 2, device = 'pdf')




tree_table = as_tibble(tree1)




# lineal tree without clades or albertii sp


ggtree(nhx_reduced,  size = 0.4, color = 'grey') +
  geom_tiplab(size = 1.5, aes(color = phylogroup, show.legend = FALSE)) +
  scale_colour_manual(values = c(
    "#E41A1C", # A
    "#377EB8", # B1
    "#4DAF4A", # B2
    "#984EA3", # C
    "#FF7F00", # D
    "#FFFF33", # E
    "#A65628", # F
    "grey90", 
    "#F781BF",  # G
    'grey'
  )) + 
  guides(colour = guide_legend(override.aes = list(size = 6, pch = 18))) 

ggsave(file = here('Rplots', 'rtree_phylogroups_noClades_lineal.pdf'), width = 100, height = 130, units = 'mm', scale = 2, device = 'pdf')



# Tree + worm imaging -----------------------------------------------------

# let's prepare the data table for the tree first

# this is for the resistance index
# there are some duplicated elements, but it seems for most of them
# there is not a big difference (for others YES, be aware)
dupstr = resistance[duplicated(resistance$Strain),] %>% 
  arrange(Strain) %>%  select(Strain) %>% t %>% as.character()

resistance %>% 
  filter(Strain %in% dupstr) %>% 
  arrange(Strain) %>% View()

# keep the first instance of the duplicated elements
res2 = resistance %>%
  select(-ID) %>% 
  rename(ID = Strain) %>% 
  distinct(ID, .keep_all = TRUE)
  

info = strain_db %>% left_join(worm.img) %>%
  # filter(!phylogroup %in% c('cladeI', 'cladeII', 'cladeIII', 'cladeIV', 'cladeV', 'E or cladeI', 'albertii')) %>% 
  drop_na(phylogroup) %>% 
  left_join(res2) %>%
  select(names, log2FC, FDR, Mean, Broadphenotype) %>%
  mutate(sig = case_when(FDR < 0.05 ~ '*',
                         FDR >= 0.05 ~ '')) %>% 
  select(-FDR) %>% 
  data.frame
  

rownames(info) = info$names
info[,1] = NULL

# df for log2FC
logFC = data.frame(info$log2FC)
rownames(logFC) = rownames(info)
colnames(logFC) = c('log2FC')

# df for sig
sig = data.frame(info$sig)
rownames(sig) = rownames(info)
colnames(sig) = c('sig')

# resistance
res3 = data.frame(info$Mean)
rownames(res3) = rownames(info)
colnames(res3) = c('Resistance')


p = ggtree(tree1, layout = 'fan', open.angle = 20, branch.length="none", size = 0.4, color = 'grey') +
  geom_tiplab(size = 1.5, aes(color = phylogroup, angle = angle, show.legend = FALSE)) +
  scale_colour_manual(values = c(
    "#E41A1C", # A
    "grey90",
    "#377EB8", # B1
    "#4DAF4A", # B2
    "#984EA3", # C
    "grey90",
    "grey90",
    "grey90",
    "grey90",
    "grey90",
    "#FF7F00", # D
    "#FFFF33", # E
    "#FFFF33", # E or clade I
    "#A65628", # F
    "grey90", 
    "#F781BF",  # G
    'grey'
  )) + guides(colour = guide_legend(override.aes = list(size = 6, pch = 18)))

# add imaging 
p1 = gheatmap(p, (logFC), offset = 4,  width=.1,
         colnames_angle = 270, hjust = -0.2, colnames_offset_y = 0.25) +
  scale_fill_viridis_c(option="B", name = 'log2FC\n imaging')


library(ggnewscale)
p2 = p1 + new_scale_fill()

p2 = gheatmap(p2, res3, offset = 8,  width=.1,
         colnames_angle = 270, hjust = -0.2, colnames_offset_y = 0.25) +
  scale_fill_viridis_c(name = 'resistance')


p3 = p2 + new_scale_fill()

gheatmap(p3, sig, offset = 7.7, width=.02,
         colnames_angle = 270, hjust = -0.2, colnames_offset_y = .25) +
  scale_fill_viridis_d(option = "D", name = "sig")

ggsave(file = here('Rplots', 'rtree_phylog_log2FC.pdf'), width = 120, height = 120, units = 'mm', scale = 2, device = 'pdf')








