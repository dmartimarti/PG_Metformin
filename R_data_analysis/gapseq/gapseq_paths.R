
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






# in depth analysis -------------------------------------------------------










