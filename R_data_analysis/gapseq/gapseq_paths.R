
# libraries ---------------------------------------------------------------



library(tidyverse)
library(readr)
library(cowplot)
library(here)
library(ComplexHeatmap)
library(circlize)


theme_set(theme_cowplot(15))



# load data ---------------------------------------------------------------



# get the genome names
gnm_names = list.files(path = ".", pattern = ".tbl") %>% 
  str_sub(start = 1, end = -18)

# get the files names
files_list = list.files(path = ".", pattern = ".tbl")


# loop over the files and get them in a huge table
genome_paths = tibble()
for (i in 1:length(gnm_names)){
  temp = read_delim(files_list[i], 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE, skip = 2) %>% 
    mutate(Genome = gnm_names[i])
  
  genome_paths = bind_rows(genome_paths, temp)
}


# load metadata 

library(readxl)
metadata = read_excel("~/Documents/MRC_postdoc/Pangenomic/metadata/MAIN_metadata.xlsx")


# how many genomes are we lacking ? 

# filter metadata for the bugs we want to analyse: AUS and ECOREF
meta_filt = metadata %>% 
  mutate(Genome = str_sub(fasta, start = 1, end = -7), .before = ID) %>% 
  filter(Discard == 'No') %>% 
  filter(Origin %in% c('AUS', 'ECOREF')) %>% 
  distinct(Genome, .keep_all = T) 


# save the missing genome list in a csv file
meta_filt %>% 
  filter(!(Genome %in% gnm_names)) %>% 
  write_csv('../missing_genomes.csv')
  





# exploration plots -------------------------------------------------------

# first filter the genome list we want to analyse

genome_paths = genome_paths %>% filter(Genome %in% meta_filt$Genome)


# calculate how many pathways per genome are in our dataset


## histograms ####

genome_paths %>% 
  filter(Completeness == 100) %>% 
  select(Name, Genome) %>% 
  group_by(Genome) %>% 
  count() %>% 
  filter(n > 400) %>% 
  ggplot(aes(n)) +
  geom_histogram(color = 'black', fill = 'grey50', bins = 30) +
  labs(
    x = 'Number of complete pathways in genome',
    y = 'Number of genomes'
  )

ggsave("../exploration/number_of_paths.pdf", height = 7, width = 9)


genome_paths %>% 
  filter(Completeness > 0) %>% 
  select(Name, Genome) %>% 
  group_by(Genome) %>% 
  count() %>% 
  filter(n > 400) %>%
  ggplot(aes(n)) +
  geom_histogram(color = 'black', fill = 'grey50', bins = 30) +
  labs(
    x = 'Number of pathways in genome',
    y = 'Number of genomes',
    caption = 'Pathways not complete (<100% of coverage) are also included in this plot'
  )
ggsave("../exploration/number_of_paths_all.pdf", height = 7, width = 9)

## complete paths matrix ####
paths_pa = genome_paths %>% 
  filter(Completeness == 100) %>% 
  select(Name, Genome) %>% 
  mutate(presence = 1) %>% 
  pivot_wider(names_from = Name, values_from = presence, values_fill = 0)

pahts_names = names(paths_pa)
gnm_names = paths_pa$Genome

paths_matrix = paths_pa %>% 
  select(-Genome) %>% 
  as.matrix() 
  
rownames(paths_matrix) = gnm_names


col_fun = colorRamp2(c(0, 1), c("white", "#0949AB"))
col_fun(seq(-3, 3))
Heatmap(paths_matrix, 
        col = col_fun,
        name = "Completeness")


quartz.save(file = '../exploration/Pathways_heatmap.pdf',
            type = 'pdf', height = 70, width = 90)


Heatmap(paths_matrix, 
        col = col_fun,
        name = "Completeness",
        show_row_names = FALSE,
        show_column_names = FALSE)

quartz.save(file = '../exploration/Pathways_heatmap_noNames.pdf',
            type = 'pdf', height = 9, width = 12)


## difference paths matrix  ####

paths_pa = genome_paths %>% 
  filter(Completeness == 100) %>% 
  select(Name, Genome) %>% 
  mutate(presence = 1) %>% 
  pivot_wider(names_from = Name, values_from = presence, values_fill = 0)

pahts_names = names(paths_pa)
gnm_names = paths_pa$Genome

paths_matrix = paths_pa %>% 
  select(-Genome) %>% 
  as.matrix() 

rownames(paths_matrix) = gnm_names

# total genomes
tot_gnms = dim(paths_matrix)[1]

# filter the paths that are NOT 100% complete in at least 99% of genomes
# (core pathway)
redux_matrix = paths_matrix[,colSums(paths_matrix) < tot_gnms*0.99]
dim(redux_matrix)


col_fun = colorRamp2(c(0, 1), c("white", "#0949AB"))
Heatmap(redux_matrix, 
        col = col_fun,
        name = "Completeness")

quartz.save(file = '../exploration/Pathways_heatmap_differential.pdf',
            type = 'pdf', height = 50, width = 70)

Heatmap(redux_matrix, 
        col = col_fun,
        name = "Completeness",
        show_row_names = FALSE,
        show_column_names = FALSE)

quartz.save(file = '../exploration/Pathways_heatmap_differential_noNames.pdf',
            type = 'pdf', height = 6, width = 9)



## all paths matrix ####

paths_pa = genome_paths %>% 
  filter(Completeness > 0) %>% 
  mutate(presence = Completeness/100) %>% 
  select(Name, Genome, presence) %>% 
  # mutate(presence = 1) %>% 
  pivot_wider(names_from = Name, values_from = presence, values_fill = 0)

pahts_names = names(paths_pa)
gnm_names = paths_pa$Genome

paths_matrix = paths_pa %>% 
  select(-Genome) %>% 
  as.matrix() 

rownames(paths_matrix) = gnm_names

Heatmap(paths_matrix, 
        col = col_fun,
        name = "Completeness")

quartz.save(file = '../exploration/ALL_Pathways_heatmap.pdf',
            type = 'pdf', height = 70, width = 190)

Heatmap(paths_matrix, 
        col = col_fun,
        name = "Completeness",
        show_row_names = FALSE,
        show_column_names = FALSE)

quartz.save(file = '../exploration/ALL_Pathways_heatmap_noNames.pdf',
            type = 'pdf', height = 9, width = 12)






# metabolites -------------------------------------------------------

# read metabolites extracted from the R objects (see the other script)
metab = read_csv("~/Documents/MRC_postdoc/Pangenomic/pangenome_analysis/ALL/phylo_analysis/gapseq/model_metabolites.csv")


gnm_valid = meta_filt %>% pull(Genome)

unique(metab$comp)

gnm_remove = c('SPC_3.1','OP50',
               'NT12226_337','NT12229','NT12169_305',
               '157')

metab %>% 
  filter(!(model %in% gnm_remove)) %>% 
  filter(model %in% gnm_valid) %>% 
  group_by(model) %>% 
  count() %>% 
  ggplot(aes(n)) +
  geom_histogram(color = 'black', fill = 'grey50', bins = 30) +
  labs(
    x = 'Molecules',
    y = 'Count'
  )

ggsave("../exploration/molecules_histogram.pdf", height = 7, width = 9)
  



metab %>% distinct(comp)

## periplasm metabs ####

# simple version with ggplot
metab %>% 
  filter(model %in% gnm_valid) %>% 
  filter(comp == 'p0') %>% 
  mutate(val = 1) %>% 
  ggplot(aes(met_name, model, fill = val)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90))

# complexHeatmap
mets_pa = metab %>% 
  filter(model %in% gnm_valid) %>% 
  filter(comp == 'p0') %>% 
  select(met_name, model) %>% 
  mutate(presence = 1) %>% 
  pivot_wider(names_from = met_name, values_from = presence, values_fill = 0)

mets_names = names(mets_pa)
gnm_names = mets_pa$model

mets_matrix = mets_pa %>% 
  select(-model) %>% 
  as.matrix() 

rownames(mets_matrix) = gnm_names


col_fun = colorRamp2(c(0, 1), c("white", "#0949AB"))

Heatmap(mets_matrix, 
        col = col_fun,
        name = "Presence")


quartz.save(file = '../exploration/Metabs_periplasm_heatmap.pdf',
            type = 'pdf', height = 80, width = 9)





## external metabs ####


metab %>% 
  filter(model %in% gnm_valid) %>% 
  filter(comp == 'e0') %>% 
  mutate(val = 1) %>% 
  ggplot(aes(met_name, model, fill = val)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90))


# complexHeatmap
mets_pa = metab %>% 
  filter(model %in% gnm_valid) %>% 
  filter(comp == 'e0') %>% 
  select(met_name, model) %>% 
  mutate(presence = 1) %>% 
  pivot_wider(names_from = met_name, values_from = presence,
              values_fn = mean, values_fill = 0)

mets_names = names(mets_pa)
gnm_names = mets_pa$model

mets_matrix = mets_pa %>% 
  select(-model) %>% 
  as.matrix() 

rownames(mets_matrix) = gnm_names


col_fun = colorRamp2(c(0, 1), c("white", "#0949AB"))

Heatmap(mets_matrix, 
        col = col_fun,
        name = "Presence")


quartz.save(file = '../exploration/Metabs_ext_heatmap.pdf',
            type = 'pdf', height = 80, width = 30)
