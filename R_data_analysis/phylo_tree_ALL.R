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
# metadata
# metadata = read_excel("D:/MRC_Postdoc/Pangenomic/pangenome_analysis/metadata/MAIN_metadata.xlsx", 
#                       sheet = "metadata")
# metadata <- read_excel("data/MAIN_metadata.xlsx")

metadata = read_excel("data/Metadata_phenotype_extended.xlsx")


# load tree
# tree = read.tree('D:/MRC_Postdoc/Pangenomic/pangenome_analysis/ALL/phylo_analysis/panaroo_results_noEVO/core_fast_tree.treefile')
tree = read.tree('data/core_fast_tree.treefile')

## WORM
# summary statistics from worm imaging

worm.img =  read_excel("data/ALL_worm_stats.xlsx", 
                       sheet = "Stats_per_strain")



## BACTERIA

resistance = read_csv("D:/MRC_Postdoc/Pangenomic/pangenome_analysis/ALL/bacterial_growth/bacterial_growth_ALL.csv")  %>% 
  select(-`...1`)


# # # # # # # # # # # # # # # #
# general tree ------------------------------------------------------------
# # # # # # # # # # # # # # # # 


# phylogroup info
sp_names = c(tree$tip.label)
df = data.frame(sp_names)
colnames(df) = 'label'


metadata = metadata %>% 
  mutate(names = str_sub(fasta, 1, -7), .before = Strainname)

# check
metadata %>% 
  filter(names %in% sp_names)

# tree %>% as_tibble() %>% full_join(df) %>% view()


# add phylogroup info and fix NT12139 strain
df = df %>% left_join(metadata %>% 
                        select(names, phylogroup) %>%
                        rename(label = names)) 
unique(df$phylogroup)


# root tree and specify groups 
grp = list(
  'A' = df %>% filter(phylogroup == 'A') %>% pull(label),
  'B1' = df %>% filter(phylogroup == 'B1') %>% pull(label),
  'B2' = df %>% filter(phylogroup == 'B2') %>% pull(label),
  'C' = df %>% filter(phylogroup == 'C') %>% pull(label),
  'G' = df %>% filter(phylogroup == 'G') %>% pull(label),
  'F' = df %>% filter(phylogroup == 'F') %>% pull(label),
  'E' = df %>% filter(phylogroup == 'E') %>% pull(label),
  'D' = df %>% filter(phylogroup == 'D') %>% pull(label),
  'E or cladeI' = df %>% filter(phylogroup == 'E or cladeI') %>% pull(label) ,
  'cladeI' = df %>% filter(phylogroup == 'cladeI') %>% pull(label)
)

tree = tree %>% groupOTU(grp, 'phylogroup')


tree_phylogroup = function(tree) {
  
  # phylogroup info
  sp_names = c(tree$tip.label)
  df = data.frame(sp_names)
  colnames(df) = 'label'
  
  # add phylogroup info and fix NT12139 strain
  df = df %>% left_join(metadata %>% 
                          select(names, phylogroup) %>%
                          rename(label = names)) 
  unique(df$phylogroup)
  
  
  # root tree and specify groups 
  grp = list(
    'A' = df %>% filter(phylogroup == 'A') %>% pull(label),
    'B1' = df %>% filter(phylogroup == 'B1') %>% pull(label),
    'B2' = df %>% filter(phylogroup == 'B2') %>% pull(label),
    'C' = df %>% filter(phylogroup == 'C') %>% pull(label),
    'G' = df %>% filter(phylogroup == 'G') %>% pull(label),
    'F' = df %>% filter(phylogroup == 'F') %>% pull(label),
    'E' = df %>% filter(phylogroup == 'E') %>% pull(label),
    'D' = df %>% filter(phylogroup == 'D') %>% pull(label),
    'E or cladeI' = df %>% filter(phylogroup == 'E or cladeI') %>% pull(label) ,
    'cladeI' = df %>% filter(phylogroup == 'cladeI') %>% pull(label)
  )
  
  tree = tree %>% groupOTU(grp, 'phylogroup')
  return(tree)
}

tree_phylogroup(tree = tree)


# 
ggtree(tree, layout = 'fan', size = 0.4)+ aes(color = phylogroup) +
  geom_tiplab(size = 1.5, aes(color = phylogroup, angle = angle, show.legend = FALSE)) +
  scale_colour_manual(values = c(
    "#E41A1C", # A
    "#377EB8", # B1
    "#4DAF4A", # B2
    "#984EA3", # C
    'grey85',    # cladeI
    "#FF7F00", # D
    "#FFFF33", # E
    "#FFFF33", # E or clade I
    "#A65628", # F
    "#F781BF"  # G  
  )) + 
  guides(colour = guide_legend(override.aes = list(size = 6, pch = 18)))


ggsave(file = here('exploration', 'rtree_phylogroups_distances.pdf'), width = 100, height = 100, units = 'mm', scale = 2, device = 'pdf')


# plot for review (lab Animal) --------------------------------------------


ggtree(tree, layout = 'fan', size = 0.4, branch.length = 'none')

ggsave(file = here('exploration', 'rtree_nolabels.eps'), width = 100, height = 100, units = 'mm', scale = 2, device = 'pdf')






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


info = worm.img %>%
  select(ID:FDR_stars) %>%  # select only relevant categories
  left_join(metadata) %>%
  filter(Discard == 'No') %>% 
  filter(!phylogroup %in% c('cladeI', 'cladeII', 'cladeIII', 'cladeIV', 
                            'cladeV', 'E or cladeI', 'albertii')) %>%
  drop_na(phylogroup) %>% 
  select(-Plate) %>% 
  left_join(res2) %>%
  select(names, log2FC, FDR, Mean = Bact_score_mean, Broadphenotype, phylogroup) %>%
  mutate(sig = case_when(FDR < 0.05 ~ '*',
                         FDR >= 0.05 ~ '')) %>% 
  select(-FDR) %>%
  mutate(Mean = case_when(Mean > 3 ~ 2.5,
                          Mean <= 20 ~ Mean)) %>% 
  mutate(Broadphenotype = replace_na(Broadphenotype, 'Unknown')) %>%
  mutate(Broadphenotype = case_when(Broadphenotype == 'Commensalstrain' ~ 'Commensal strain',
                                    TRUE ~ Broadphenotype)) %>% 
  distinct(names, .keep_all = T) %>% 
  data.frame


rownames(info) = info$names
info[,1] = NULL

unique(info$Broadphenotype)


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


p = ggtree(tree, layout = 'fan', open.angle = 25, branch.length="none", size = 0.5, color = 'grey') +
  geom_tiplab(size = 1.5, aes(color = phylogroup, angle = angle, show.legend = FALSE)) +
  scale_colour_manual(values = c(
    "#E41A1C", # A
    "#377EB8", # B1
    "#4DAF4A", # B2
    "#984EA3", # C
    'grey85',    # cladeI
    "#FF7F00", # D
    "#FFFF33", # E
    "#FFFF33", # E or clade I
    "#A65628", # F
    "#F781BF"  # G  
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

ggsave(file = here('exploration', 'rtree_phylog_log2FC.pdf'), 
       width = 180, height = 180, units = 'mm', scale = 2, device = 'pdf')






# # # # # # # # # # # # # # # # # # # # # # 
# Tree + worm 2 imaging -----------------------------------------------------

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


# calculate distances from metadata
worm.img.dist = worm.img %>%  
  select(ID:FDR_stars) %>%  # select only relevant categories
  left_join(metadata) %>%
  filter(Discard == 'No') %>% 
  mutate(Annotation_0mM = case_when(Annotation_0mM == 'normal' | Annotation_0mM == 'NB' ~ 0,
                                    Annotation_0mM == 'biofilm' ~ 1,
                                    Annotation_0mM == 'super_bio' | Annotation_0mM == 'super_biofilm' ~ 2),
         Annotation_50mM = case_when(Annotation_50mM == 'normal' ~ 0,
                                     Annotation_50mM == 'biofilm' ~ 1,
                                     Annotation_50mM == 'super_bio' |Annotation_50mM == 'super_biofilm' ~ 2),
         distance = Annotation_50mM - Annotation_0mM) %>% 
  drop_na(ID) %>% 
  arrange(desc(distance)) %>% 
  select(ID, Broadphenotype, phylogroup, distance)

# check that all distances are ok
is.na(worm.img.dist$distance)
worm.img.dist$distance

info = worm.img %>%
  select(ID:FDR_stars) %>%  # select only relevant categories
  left_join(metadata) %>%
  filter(Discard == 'No') %>% 
  left_join(worm.img.dist %>% 
              select(ID, distance)) %>% 
  filter(!phylogroup %in% c('cladeI', 'cladeII', 'cladeIII', 'cladeIV', 'cladeV', 'E or cladeI', 'albertii')) %>%
  drop_na(phylogroup) %>% 
  select(-Plate) %>% 
  left_join(res2) %>%
  select(names, log2FC, FDR, distance, Mean = Bact_score_mean, Broadphenotype) %>%
  mutate(sig = case_when(FDR < 0.05 ~ '*',
                         FDR >= 0.05 ~ ''),
         distance = factor(distance, levels = c(-2,-1,0,1,2))) %>% 
  select(-FDR) %>%
  mutate(Mean = case_when(Mean > 3 ~ 2.5,
                          Mean <= 20 ~ Mean)) %>% 
  mutate(Broadphenotype = replace_na(Broadphenotype, 'Unknown')) %>%
  mutate(Broadphenotype = case_when(Broadphenotype == 'Commensalstrain' ~ 'Commensal strain',
                                    TRUE ~ Broadphenotype)) %>% 
  distinct(names, .keep_all = T) %>% 
  data.frame




rownames(info) = info$names
info[,1] = NULL

# df for log2FC
logFC = data.frame(info$log2FC)
rownames(logFC) = rownames(info)
colnames(logFC) = c('acs2_FC')

# df for sig
sig = data.frame(info$sig)
rownames(sig) = rownames(info)
colnames(sig) = c('sig')

# resistance
res3 = data.frame(info$Mean)
rownames(res3) = rownames(info)
colnames(res3) = c('Bacterial_score')

# phenotype
pheno = data.frame(info$Broadphenotype)
rownames(pheno) = rownames(info)
colnames(pheno) = 'Phenotype'

# distance
distance = data.frame(info$distance)
rownames(distance) = rownames(info)
colnames(distance) = 'Distance'




p = ggtree(tree, layout = 'fan', open.angle = 25, branch.length="none", size = 0.5, color = 'grey') +
  geom_tiplab(size = 1.5, aes(color = phylogroup, angle = angle, show.legend = FALSE)) +
  scale_colour_manual(values = c(
    "#E41A1C", # A
    "#377EB8", # B1
    "#4DAF4A", # B2
    "#984EA3", # C
    'grey85',    # cladeI
    "#FF7F00", # D
    "#FFFF33", # E
    "#FFFF33", # E or clade I
    "#A65628", # F
    "#F781BF"  # G  
  )) + guides(colour = guide_legend(override.aes = list(size = 6, pch = 18)))

# add imaging 
p1 = gheatmap(p, (logFC), offset = 4,  width=.1,
              colnames_angle = 270, hjust = -0.2, colnames_offset_y = 0.25) +
  scale_fill_viridis_c(option="B", name = 'acs-2 FC')


library(ggnewscale)
p2 = p1 + new_scale_fill()

p2 = gheatmap(p2, res3, offset = 8,  width=.1,
              colnames_angle = 270, hjust = -0.2, colnames_offset_y = 0.25) +
  scale_fill_viridis_c(name = 'Bacterial score')


p3 = p2 + new_scale_fill()

p3 = gheatmap(p3, sig, offset = 8.5, width=.02,
              colnames_angle = 270, hjust = -0.2, colnames_offset_y = .25) +
  scale_fill_manual(values = c('black','grey'), name = "sig")


p4 = p3 + new_scale_fill()

p4 = gheatmap(p4, pheno, offset = 11.1, width=0.06,
              colnames_angle = 270, hjust = -0.2, colnames_offset_y = .25) +
  scale_fill_manual(values = c('#27C73C',
                               '#C73F16',
                               '#0228C7',
                               'grey50'), name = "Phenotype")

p5 = p4 + new_scale_fill()

gheatmap(p5, distance, offset = 13,  width=.1,
         colnames_angle = 270, hjust = -0.2, colnames_offset_y = 0.25) +
  scale_fill_manual(name = 'Distance',
                    values = c('#FB8A21',
                               '#E31EE0',
                               '#2D82FA',
                               '#1DA82F',
                               '#FFD424'))

ggsave(file = here('exploration', 'rtree_phylog_biofilm.pdf'), width = 180, height = 180, units = 'mm', scale = 2, device = 'pdf')





# tree reordered --------------------



# df for log2FC
logFC = data.frame(info$log2FC)
rownames(logFC) = rownames(info)
colnames(logFC) = c('Pacs2::GFP FC')

# df for sig
sig = data.frame(info$sig)
rownames(sig) = rownames(info)
colnames(sig) = c('Significance')

# resistance
res3 = data.frame(info$Mean)
rownames(res3) = rownames(info)
colnames(res3) = c('Bacterial Score')

# phenotype
pheno = data.frame(info$Broadphenotype)
rownames(pheno) = rownames(info)
colnames(pheno) = 'Bacterial Phenotype'

# distance
distance = data.frame(info$distance)
rownames(distance) = rownames(info)
colnames(distance) = 'Phenotype distance'



p = ggtree(tree, layout = 'fan', open.angle = 25, branch.length="none", size = 0.5, color = 'grey') +
  geom_tiplab(size = 1.5, aes(color = phylogroup, angle = angle, show.legend = FALSE)) +
  scale_colour_manual(values = c(
    "#E41A1C", # A
    "#377EB8", # B1
    "#4DAF4A", # B2
    "#984EA3", # C
    'grey85',    # cladeI
    "#FF7F00", # D
    "#FFFF33", # E
    "#FFFF33", # E or clade I
    "#A65628", # F
    "#F781BF"  # G  
  )) + guides(colour = guide_legend(override.aes = list(size = 6, pch = 18)))



# 1. bacterial phenotype
p1 = gheatmap(p, pheno, offset = 4, width=0.1,
              colnames_angle = 0, hjust = -0.2, colnames_offset_y = .25) +
  scale_fill_manual(values = c('#27C73C',
                               '#C73F16',
                               '#0228C7',
                               'grey50'), name = "Bacterial Phenotype")

# 2. bacterial resistance score
p2 = p1 + new_scale_fill()

p2 = gheatmap(p2, res3, offset = 7.2,  width=.1,
              colnames_angle = 0, hjust = -0.2, colnames_offset_y = 0.25) +
  scale_fill_viridis_c(name = 'Bacterial score')


# 3. distance
p3 = p2 + new_scale_fill()

p3 = gheatmap(p3, distance, offset = 10.4,  width=.1,
              colnames_angle = 0, hjust = -0.2, colnames_offset_y = 0.25) +
  scale_fill_manual(name = 'Phenotype distance',
                    values = c( '#B3FA41', # very green
                                '#43E6BD', # algae green
                                '#5963FA', # blue
                                '#FF8D38', # orange 
                                '#E64040'  # red
                    ) 
  )

# 4. Pacs-2::GFP imaging 

p4 = p3 + new_scale_fill()

p4 = gheatmap(p4, logFC, offset = 13.1,  width=.1,
              colnames_angle = 0, hjust = -0.2, colnames_offset_y = 0.25) +
  scale_fill_viridis_c(option="B", name = 'Pacs2::GFP FC')


# 5. significance for FC 

p5 = p4 + new_scale_fill()

p5 = gheatmap(p5, sig, offset = 17.6, width=.02,
              colnames_angle = 0, hjust = -0.2, colnames_offset_y = .25) +
  scale_fill_manual(values = c('black','grey'), name = "Significance")

rotate_tree(p5, 90)


ggsave(file = here('exploration', 'rtree_phylog_biofilm.pdf'), width = 180, height = 180, units = 'mm', scale = 2, device = 'pdf')




# Filipe's grant ----------------------------------------------------------

# info = worm.img %>%
#   select(ID:FDR_stars) %>%  # select only relevant categories
#   left_join(metadata) %>%
#   filter(Discard == 'No') %>%
#   left_join(worm.img.dist %>%
#               select(ID, distance)) %>%
#   filter(!phylogroup %in% c('cladeI', 'cladeII', 'cladeIII', 'cladeIV', 'cladeV', 'E or cladeI', 'albertii')) %>%
#   drop_na(phylogroup) %>%
#   select(-Plate) %>%
#   left_join(res2) %>%
#   select(names, log2FC, FDR, distance, Mean = Bact_score_mean, Broadphenotype) %>%
#   mutate(sig = case_when(FDR < 0.05 ~ '*',
#                          FDR >= 0.05 ~ ''),
#          distance = factor(distance, levels = c(-2,-1,0,1,2))) %>%
#   select(-FDR) %>%
#   mutate(Mean = case_when(Mean > 3 ~ 2.5,
#                           Mean <= 20 ~ Mean)) %>%
#   mutate(Broadphenotype = replace_na(Broadphenotype, 'Unknown')) %>%
#   mutate(Broadphenotype = case_when(Broadphenotype == 'Commensalstrain' ~ 'Commensal strain',
#                                     TRUE ~ Broadphenotype)) %>%
#   distinct(names, .keep_all = T) %>%
#   data.frame


info = worm.img %>%
  select(ID:FDR_stars) %>%  # select only relevant categories
  left_join(metadata) %>%
  filter(Discard == 'No') %>%
  left_join(worm.img.dist %>%
              select(ID, distance)) %>%
  filter(!phylogroup %in% c('cladeI', 'cladeII', 'cladeIII', 'cladeIV', 'cladeV', 'E or cladeI', 'albertii')) %>%
  drop_na(phylogroup) %>%
  select(-Plate) %>%
  left_join(res2 %>% select(ID, Bact_score_mean)) %>% 
  mutate(Bact_score_mean = case_when(Bact_score_mean > 3 ~ 3,
                                     TRUE ~ Bact_score_mean)) %>% 
  mutate(extended_phenotype = case_when(is.na(extended_phenotype) ~ 'Unknown',
                                        TRUE ~ extended_phenotype)) %>% 
  # mutate(Broadphenotype = replace_na(Broadphenotype, 'Unknown')) %>%
  # mutate(Broadphenotype = case_when(Broadphenotype == 'Commensalstrain' ~ 'Commensal strain',
  #                                   TRUE ~ Broadphenotype)) %>% 
  distinct(names, .keep_all = T) %>% 
  data.frame


rownames(info) = info$names
info[,1] = NULL


# df for log2FC
logFC = data.frame(info$log2FC)
rownames(logFC) = rownames(info)
colnames(logFC) = c('Microbe-Drug-Host')

# df for sig
base = data.frame(info$Met_0mM)
rownames(base) = rownames(info)
colnames(base) = c('Microbe-Host')

# resistance
res3 = data.frame(info$Bact_score_mean)
rownames(res3) = rownames(info)
colnames(res3) = c('Microbial Growth')

# phenotype
pheno = data.frame(info$extended_phenotype)
rownames(pheno) = rownames(info)
colnames(pheno) = 'Microbial Ecological Niche'

# # distance
# distance = data.frame(info$distance)
# rownames(distance) = rownames(info)
# colnames(distance) = 'Phenotype distance'

tree$tip.label[!(tree$tip.label %in% row.names(logFC))]

to_drop = 
  tree$tip.label[!(tree$tip.label %in% row.names(logFC))]

tree_reduced = drop.tip(tree, to_drop)

tree_reduced = tree_phylogroup(tree = tree_reduced)

p = ggtree(tree_reduced, layout = 'fan', open.angle = 25, branch.length="none", size = 0.5, color = 'grey') +
  geom_tiplab(size = 1.3, aes(color = phylogroup, angle = angle, show.legend = FALSE)) +
  scale_colour_manual(values = c(
    "#E41A1C", # A
    "#377EB8", # B1
    "#4DAF4A", # B2
    "#984EA3", # C
    # 'grey85',    # cladeI
    "#FF7F00", # D
    "#FFFF33", # E
    # "#FFFF33", # E or clade I
    "#A65628", # F
    "#F781BF"  # G  
  )) + guides(colour = guide_legend(override.aes = list(size = 6, pch = 18)))





# 1. bacterial phenotype
p1 = gheatmap(p, pheno, offset = 4, width=0.1,
              colnames_angle = 0, hjust = -0.06, colnames_offset_x = 0.25, color=NA) +
  scale_fill_manual(values = c('#1D6EDB',
                               '#1CDB06',
                               '#DB1D97',
                               '#DB9612',
                               'grey50'), name = "Bacterial Phenotype")

# 2. bacterial resistance score
p2 = p1 + new_scale_fill()

p2 = gheatmap(p2, res3, offset = 7.2,  width=.1,
              colnames_angle = 0, hjust = -0.2, colnames_offset_y = 0.25,color=NA) +
  scale_fill_viridis_c(name = 'Bacterial score', direction = -1)


# 3. Fold change of worm imaging
p3 = p2 + new_scale_fill()

p3 = gheatmap(p3, logFC, offset = 13.1,  width=.1,
              colnames_angle = 0, hjust = -0.2, colnames_offset_y = 0.25,color=NA) +
  scale_fill_viridis_c(option="D",direction = -1, name = 'Pacs2::GFP log2 Fold Change')

# 4. GFP control values
p4 = p3 + new_scale_fill()

p4 = gheatmap(p4, base, offset = 10.1,  width=.1,
              colnames_angle = 0, hjust = -0.2, colnames_offset_y = 0.25,color=NA) +
  scale_fill_viridis_c(option="B", direction = -1, name = 'Pacs2::GFP Control')


rotate_tree(p4, 90)


ggsave(file = here('exploration', 'grant_tree.pdf'),
       width = 200, height = 200, units = 'mm', scale = 2, 
       device = 'pdf')




# iTOL files --------------------------------------------------------------


# let's create an iTOL annotation file to colour the labels in the tree

library(readxl)
metadata = read_excel("tree_files/MAIN_metadata.xlsx")



### IMPORTANT

strains = read_delim("tree_files/strains.tsv", 
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE)



strains = strains %>% 
  mutate(extended_phenotype = 
           case_when(Phenotype == 'Laboratory strain' ~ 'Laboratory strain',
                     str_detect(`Isolation origin`, 'Human') ~ 'Human commensal',
                     str_detect(`Isolation origin`, 'human') ~ 'Human commensal',
                     Phenotype == 'Commensal strain' ~ 'Animal commensal',
                     Phenotype == 'Pathogenic strain' ~ 'Pathogenic strain',
                     is.na(Phenotype) ~ 'Unknown',
                     Phenotype == 'Unknown' ~ 'Unknown',
                     Phenotype == 'Evolution experiment' ~ 'Evolution experiment'),
         .before = `Parental Strain`) 


metadata %>% 
  left_join(
    strains %>% 
      select(ID = `Strain Identifier`, extended_phenotype)
  ) %>% 
  mutate(extended_phenotype = 
           case_when(is.na(extended_phenotype) ~ Broadphenotype,
                     TRUE ~ extended_phenotype)) %>% 
  mutate(extended_phenotype = str_replace_all(extended_phenotype, 
                                              'Commensal strain',
                                              'Human commensal'),
         extended_phenotype = case_when(
           Broadphenotype == 'Pathogenic strain' ~ 'Pathogenic strain',
           TRUE ~ extended_phenotype)) %>% 
  select(fasta:Broadphenotype, extended_phenotype, everything()) %>% 
  write.xlsx('Metadata_phenotype_extended.xlsx')

metadata = read_excel('Metadata_phenotype_extended.xlsx')

# "#E41A1C", # A
# "#377EB8", # B1
# "#4DAF4A", # B2
# "#984EA3", # C
# 'E8E8E8',    # cladeI
# "#FF7F00", # D
# "#FFFF33", # E
# "#FFFF33", # E or clade I
# "#A65628", # F
# "#F781BF"  # G  

metadata %>% 
  filter(Discard == 'No') %>% 
  mutate(ID = str_replace_all(fasta, '.fasta', '')) %>% 
  select(ID, phylogroup) %>% 
  mutate(color = 
           case_when(phylogroup == 'A' ~ '#E41A1C',
                     phylogroup == 'B1' ~ '#377EB8',
                     phylogroup == 'B2' ~ '#4DAF4A',
                     phylogroup == 'C' ~ '#984EA3',
                     phylogroup == 'cladeI' ~ '#E8E8E8',
                     phylogroup == 'D' ~ '#FF7F00',
                     phylogroup == 'E' ~ '#FFFF33',
                     phylogroup == 'E or clade I' ~ '#FFFF33',
                     phylogroup == 'F' ~ '#A65628',
                     phylogroup == 'G' ~ '#F781BF',
           )
  ) %>% 
  mutate(type = 'label', 
         font = 'bold', 
         size = 1) %>% 
  select(ID, type, color, font, size) %>% 
  write_delim('label_colours.txt')



metadata %>% 
  filter(Discard == 'No') %>% 
  mutate(ID = str_replace_all(fasta, '.fasta', '')) %>% 
  select(ID, extended_phenotype) %>% 
  drop_na() %>% 
  mutate(color = 
           case_when(
             extended_phenotype == 'Human commensal' ~ '#219FEB',
             extended_phenotype == 'Animal commensal' ~ '#5CEB09',
             extended_phenotype == 'Pathogenic strain' ~ '#EB9015',
             extended_phenotype == 'Laboratory strain' ~ '#E821EB'
           )) %>% 
  select(ID, color, extended_phenotype) %>% 
  write_delim('colorstrip_vars.txt', col_names = FALSE, quote = 'none')
  



