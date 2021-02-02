
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


theme_set(theme_light())

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
# adsaresistance = read_excel("200601_2v3_top50_sensitivePG.xlsx", 
                                            # sheet = "111")
resistance = read_excel("D:/MRC_Postdoc/Pangenomic/biolog/PG_growth_biospa/Growth_resistance_summaryStats.xlsx", 
                                             sheet = "unweighted") %>% select(-`...1`)


# extended info about strains
strains = read_delim("strains.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)


strains %>% 
  filter(`Strain Identifier` %in% strain_db2$ID)



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
grp = list(
'A' = df %>% filter(phylogroup == 'A') %>% select(label) %>% t %>% as.character,
'albertii' = df %>% filter(phylogroup == 'albertii') %>% select(label) %>% t %>% as.character,
'B1' = df %>% filter(phylogroup == 'B1') %>% select(label) %>% t %>% as.character,
'B2' = df %>% filter(phylogroup == 'B2') %>% select(label) %>% t %>% as.character,
'C' = df %>% filter(phylogroup == 'C') %>% select(label) %>% t %>% as.character,
'G' = df %>% filter(phylogroup == 'G') %>% select(label) %>% t %>% as.character,
'F' = df %>% filter(phylogroup == 'F') %>% select(label) %>% t %>% as.character,
'E' = df %>% filter(phylogroup == 'E') %>% select(label) %>% t %>% as.character,
'D' = df %>% filter(phylogroup == 'D') %>% select(label) %>% t %>% as.character,
'E or cladeI' = df %>% filter(phylogroup == 'E or cladeI') %>% select(label) %>% t %>% as.character,
'cladeII' = df %>% filter(phylogroup == 'cladeII') %>% select(label) %>% t %>% as.character,
'cladeV' = df %>% filter(phylogroup == 'cladeV') %>% select(label) %>% t %>% as.character,
'cladeIV' = df %>% filter(phylogroup == 'cladeIV') %>% select(label) %>% t %>% as.character,
'cladeIII' = df %>% filter(phylogroup == 'cladeIII') %>% select(label) %>% t %>% as.character,
'fergusonii' = df %>% filter(phylogroup == 'fergusonii') %>% select(label) %>% t %>% as.character,
'cladeI' = df %>% filter(phylogroup == 'cladeI') %>% select(label) %>% t %>% as.character
)

tree = tree %>% root('NT12461_14') %>% groupOTU(grp, 'phylogroup')

ggtree(tree, layout = 'fan', size = 0.4, branch.length = 'none')+ aes(color = phylogroup) +
  geom_tiplab(size = 1.5, aes(color = phylogroup, angle = angle, show.legend = FALSE)) +
  scale_colour_manual(values = c(
    "grey90",
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
    "grey90", # E or clade I
    "#A65628", # F
    "grey90", 
    "#F781BF",  # G
    'grey'
  )) + 
  guides(colour = guide_legend(override.aes = list(size = 6, pch = 18)))


ggsave(file = here('Rplots', 'rtree_phylogroups.pdf'), width = 100, height = 100, units = 'mm', scale = 2, device = 'pdf')


## remove tips from clades

to_drop = strain_db %>% 
  filter(phylogroup %in% c('cladeI', 'cladeII', 'cladeIII', 'cladeIV', 'cladeV', 'E or cladeI', 'albertii')) %>% 
  select(names) %>% t %>% as.character

# remove list elements from discarded groups
grp_reduced = grp[names(grp) %in% c('cladeI', 'cladeII', 'cladeIII', 'cladeIV', 'cladeV', 'E or cladeI', 'albertii') == FALSE]

nhx_reduced = treeio::drop.tip(tree, to_drop) %>% 
  groupOTU(grp_reduced, 'phylogroup')

ggtree(nhx_reduced, layout = 'circular', size = 0.4, branch.length = 'none') +
  aes(color = phylogroup) +
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



# lineal tree without clades or albertii sp


ggtree(nhx_reduced,  size = 0.4, aes(color = phylogroup)) +
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
# dupstr = resistance[duplicated(resistance$Strain),] %>% 
#   arrange(Strain) %>%  select(Strain) %>% t %>% as.character()
# 
# resistance %>% 
#   filter(Strain %in% dupstr) %>% 
#   arrange(Strain) %>% View()

# keep the first instance of the duplicated elements
res2 = resistance %>%
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
  mutate(Mean = case_when(Mean > 20 ~ 2.5,
                          Mean <= 20 ~ Mean)) %>% 
  mutate(Broadphenotype = replace_na(Broadphenotype, 'Unknown')) %>%  
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

# phenotype
pheno = data.frame(info$Broadphenotype)
rownames(pheno) = rownames(info)
colnames(pheno) = 'Phenotype'


p = ggtree(nhx_reduced, layout = 'fan', open.angle = 25, branch.length="none", size = 0.5, color = 'grey') +
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

p3 = gheatmap(p3, sig, offset = 7.8, width=.02,
         colnames_angle = 270, hjust = -0.2, colnames_offset_y = .25) +
  scale_fill_viridis_d(option = "D", name = "sig")


p4 = p3 + new_scale_fill()

gheatmap(p4, pheno, offset = 11.1, width=0.06,
         colnames_angle = 270, hjust = -0.2, colnames_offset_y = .25) +
  scale_fill_manual(values = c('#27C73C',
                               '#C73F16',
                               '#0228C7',
                               'grey50'), name = "Phenotype")

ggsave(file = here('Rplots', 'rtree_phylog_log2FC.pdf'), width = 180, height = 180, units = 'mm', scale = 2, device = 'pdf')














# mean vs log2FC plot with lineal models
info %>% mutate(Broadphenotype = factor(Broadphenotype)) %>% 
ggplot(aes(Mean, log2FC)) +
  geom_point() +
  geom_smooth(method= 'lm') +
  facet_wrap(~Broadphenotype)

summary(lm(Mean ~ log2FC, data = info %>% filter(Broadphenotype == 'Laboratorystrain')))



ggsave(file = here('Rplots', 'Mean_log2FC.pdf'), width = 140, height = 90, units = 'mm', scale = 2, device = 'pdf')



# Tree with categories ----------------------------------------------------

# divide Mean in quantiles and keep the 10 and 90 as the limits
lim = quantile(info$Mean, na.rm = TRUE, probs = seq(0,1,0.1))[c(3,9)]


info2 = info %>% mutate(Mean = case_when(Mean < lim[1] ~ 'Sensitive',
                                 Mean >= lim[1] & Mean <= lim[2] ~ 'Neutral',
                                 Mean > lim[2] ~ 'Resistant'))

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
res3 = data.frame(info2$Mean)
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
  scale_fill_viridis_d(option="B",name = 'resistance')


p3 = p2 + new_scale_fill()

gheatmap(p3, sig, offset = 7.8, width=.02,
         colnames_angle = 270, hjust = -0.2, colnames_offset_y = .25) +
  scale_fill_viridis_d(option = "D", name = "sig")

ggsave(file = here('Rplots', 'rtree_phylog_log2FC_categories.pdf'), width = 120, height = 120, units = 'mm', scale = 2, device = 'pdf')






# PanViz ------------------------------------------------------------------



library(PanVizGenerator)

getGO()

csvFile = read_csv("gene_presence_absence.csv")
outputDir <- here()

panviz(csvFile, location = outputDir)


