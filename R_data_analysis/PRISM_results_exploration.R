library(tidyverse)
library(tidymodels)
library(readr)
library(broom)  # devtools::install_github("tidymodels/broom")
library(cowplot)
library(ggrepel)
library(here)

# run this once

# summary_smiles = read_csv("summary_smiles.csv")
# 
# summary_smiles = summary_smiles %>% 
#   mutate(genome = as.factor(genome),
#          cluster = as.factor(cluster))

summary = read_csv("summary.csv")


summary %>% 
  count(type) %>% 
  ggplot(aes(y = n, x = fct_reorder(type, n))) +
  geom_bar(stat = 'identity') +
  theme_cowplot(12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  labs(x = NULL,
       y = 'Number of elements')

ggsave(here('exploration', 'BGCs_total.pdf'), height = 10, width = 11)


## PCA with all variants of products

pca_fit = summary_smiles %>% 
  select(-mol_weight) %>% 
  select(where(is.numeric)) %>% 
  prcomp(scale = F)

pca_fit %>%
  augment(summary_smiles) %>% # add original dataset back in
  ggplot(aes(.fittedPC1, .fittedPC2, color = genome)) + 
  geom_point(size = 1.5) +
  theme_half_open(12) + background_grid()


ggsave(here('exploration', 'PCA_compounds_TOTAL.pdf'), height = 15, width = 18)



## PCA with a redux dataset (only one product per genome per cluster)

smiles_redux = summary_smiles %>%
  distinct(genome, cluster, .keep_all = TRUE)

pca_fit = smiles_redux %>% 
  select(-mol_weight) %>% 
  select(where(is.numeric)) %>% 
  prcomp(scale = F)

pca_fit %>%
  augment(smiles_redux) %>% # add original dataset back in
  ggplot(aes(.fittedPC1, .fittedPC2)) + 
  geom_point(size = 1.5) +
  geom_text_repel(aes(label = genome), max.overlaps = 200) +
  theme_half_open(12) + 
  background_grid()

ggsave(here('exploration', 'PCA_compounds_subset.pdf'), height = 10, width = 12)


smiles_redux %>% 
  ggplot(aes(x = mol_weight)) +
  geom_histogram() +
  theme_half_open(15)

ggsave(here('exploration', 'mol_weights_dist.pdf'), height = 8, width = 9)



