library(tidyverse) # master library to deal with data frames
library(readxl) # read xlsx or xls files
library(FactoMineR) # for PCA
library(factoextra) # for PCA
library(here) # usefull to save plots in folders given a root
library(viridis) # color palette package
library(ComplexHeatmap) # yeah, complex heatmaps


strain_db = read_excel("strain_db.xlsx")

# summarise phylogroups 
phylogr = strain_db %>%
  filter(Broadphenotype != 'Evolutionexperiment') %>%
  drop_na(Assembly) %>% 
  mutate(
    phylo_corrected = ifelse(phylogroup == 'Unknown', mash_group, phylogroup),
    phylo_corrected = factor(phylo_corrected, levels = c('A', 'B2', 'B1', 'C', 'D','E', 'F','G',
                                                    'cladeI', 'cladeII', 'cladeIII', 'cladeIV','cladeV',
                                                    'E or cladeI', 'fergusonii', 'albertii', 'Non Escherichia'
                                                    ))) %>%
  mutate(phylo_corrected = fct_explicit_na(phylo_corrected)) %>%
  group_by(phylo_corrected) %>%
  summarise(N = n())

# plot phylogroups
phylogr %>%
  mutate(Total = sum(N),
         Fraction = round((N/Total)*100,1),
         y_label_pos = Fraction + 3) %>%
  ggplot(aes(x = phylo_corrected, y = Fraction, fill = phylo_corrected)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(y = y_label_pos, label = Fraction), vjust = 1.6, size = 3.5) +
  labs(x = 'phylogroups', y = '% of total') +
  theme_light() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1))


dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'phylogrous_dist_barplot.pdf'),
             height = 8, width = 8, useDingbats = FALSE)



phylogr = strain_db %>%
  filter(Broadphenotype != 'Evolutionexperiment') %>%
  drop_na(Assembly) %>% 
  mutate(
    phylo_corrected = ifelse(phylogroup == 'Unknown', mash_group, phylogroup),
    phylo_corrected = factor(phylo_corrected, levels = c('A', 'B2', 'B1', 'C', 'D','E', 'F','G',
                                                         'cladeI', 'cladeII', 'cladeIII', 'cladeIV','cladeV',
                                                         'E or cladeI', 'fergusonii', 'albertii', 'Non Escherichia'
    ))) %>%
  mutate(phylo_corrected = fct_explicit_na(phylo_corrected)) %>%
  group_by(phylo_corrected, Broadphenotype) %>%
  summarise(N = n()) %>% ungroup

phylogr %>%
  mutate(Total = sum(N),
         Fraction = round((N/Total)*100,1),
         y_label_pos = Fraction + 3) %>%
  ggplot(aes(x = phylo_corrected, y = Fraction, fill = phylo_corrected)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(y = y_label_pos, label = Fraction), vjust = 1.6, size = 3.5) +
  labs(x = 'phylogroups', y = '% of total') +
  facet_wrap(~Broadphenotype) +
  theme_light() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1))

dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'phylogrous_phenotype_dist_barplot.pdf'),
             height = 10, width = 14, useDingbats = FALSE)

