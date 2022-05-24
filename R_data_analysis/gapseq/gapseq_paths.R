
# libraries ---------------------------------------------------------------



library(tidyverse)
library(readr)
library(cowplot)
library(here)
library(ComplexHeatmap)


theme_set(theme_cowplot(15))


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



# exploration plots -------------------------------------------------------


# calculate how many pathways per genome are in our dataset

genome_paths %>% 
  filter(Completeness == 100) %>% 
  select(Name, Genome) %>% 
  group_by(Genome) %>% 
  count() %>% 
  filter(n > 400) %>% 
  ggplot(aes(n)) +
  geom_histogram(color = 'black', fill = 'grey50', bins = 30) +
  labs(
    x = 'Number of pathways in genome',
    y = 'Number of genomes'
  )

ggsave("../exploration/number_of_paths.pdf", height = 7, width = 9)


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


library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
col_fun(seq(-3, 3))
Heatmap(paths_matrix, col = col_fun)


quartz.save(file = '../exploration/Pathways_heatmap.pdf',
            type = 'pdf', height = 70, width = 90)

