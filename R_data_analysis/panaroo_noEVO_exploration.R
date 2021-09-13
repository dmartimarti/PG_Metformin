
library(tidyverse) # master library to deal with data frames
library(readxl) # read xlsx or xls files
# library(FactoMineR) # for PCA
# library(factoextra) # for PCA
library(here) # usefull to save plots in folders given a root
library(viridis) # color palette package
library(ComplexHeatmap) # yeah, complex heatmaps
library(openxlsx)
library(ggrepel)
library(patchwork)
library(cowplot)

# theme_set(theme_light())


metadata = read_excel("D:/MRC_Postdoc/Pangenomic/pangenome_analysis/metadata/MAIN_metadata.xlsx", 
                      sheet = "metadata")


# phylogroups -------------------------------------------------------------


# summarise phylogroups 
phylogr = metadata %>%
  filter(Discard == 'No') %>% 
  # filter(Broadphenotype != 'Evolutionexperiment') %>%
  drop_na(Assembly) %>% 
  mutate(
    phylo_corrected = ifelse(phylogroup == 'Unknown', mash_group, phylogroup),
    phylo_corrected = factor(phylo_corrected, levels = c('A', 'B2', 'B1', 'C', 'D','E', 'F','G', 'E or cladeI',
                                                         'cladeI', 'cladeII', 'cladeIII', 'cladeIV','cladeV',
                                                          'fergusonii', 'albertii', 'Non Escherichia'
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
  labs(x = 'Phylogroups', y = '% of total') +
  theme_cowplot(14) +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1))


dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'phylogrous_dist_barplot.pdf'),
             height = 8, width = 8, useDingbats = FALSE)





phylogr = metadata %>%
  filter(Discard == 'No') %>%
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
  facet_wrap(~Broadphenotype, ncol = 2) +
  theme_cowplot(14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'phylogrous_phenotype_dist_barplot.pdf'),
             height = 12, width = 12, useDingbats = FALSE)




# pyseer WORM ------------------------------------------------------------------


### read data ####

# ALL
all_nobio = read_delim("pyseer_output/results/worm_phenotype_ALL_no_biofilm.tsv", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)


all_worm = read_delim("pyseer_output/results/worm_phenotype_ALL.tsv", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)

# AUS strains
aus_worm = read_delim("pyseer_output/results/worm_phenotype_AUS.tsv", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)

aus_nobio = read_delim("pyseer_output/results/worm_phenotype_AUS_no_biofilm.tsv", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)

# ECOREF strains 
ecoref_nobio = read_delim("pyseer_output/results/worm_phenotype_ECOREF_no_biofilm.tsv", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)

ecoref_worm = read_delim("pyseer_output/results/worm_phenotype_ECOREF.tsv", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)



all_nobio %>% 
  arrange(variant) %>% 
  mutate(log_lrt = -log10(`lrt-pvalue`),
         log_pval = -log10(`filter-pvalue`)) %>% 
  arrange(desc(af)) %>%
  mutate(variant_single = str_split(variant, pattern = '~~~') %>% 
           map_chr(., 1)) %>% 
  mutate(
    variant_2 = case_when(
      grepl("group_",variant_single) ~ 'Unknown',
      TRUE ~ variant_single),
    labels = case_when(log_lrt > 3 ~ variant_2),
    variant = factor(variant, levels = variant)) %>% 
  mutate(labels = case_when(labels == 'Unknown' ~ '',
                           TRUE ~ labels))


# helper function to draw the Manhattan plot

manhplot.lrt = function(data, limit = 3, unk = TRUE) {
  
  # Draws a manhattan plot from a pyseer output
  # limit marks the threshold by which the funciton will add labels
  # unk: wether you want the Unknowns to be plotted or not

  data2plot = data %>% 
    arrange(variant) %>% 
    mutate(log_lrt = -log10(`lrt-pvalue`),
           log_pval = -log10(`filter-pvalue`)) %>% 
    arrange(desc(af)) %>%
    mutate(variant_single = str_split(variant, pattern = '~~~') %>% 
             map_chr(., 1)) %>% 
    mutate(
      variant_2 = case_when(
        grepl("group_",variant_single) ~ 'Unknown',
        TRUE ~ variant_single),
      labels = case_when(log_lrt > limit ~ variant_2),
      variant = factor(variant, levels = variant)) 
  
  if (unk == FALSE){
    data2plot = data2plot %>% 
      mutate(labels = case_when(labels == 'Unknown' ~ '',
                                TRUE ~ labels))
  } else if (unk == TRUE) {
    data2plot = data2plot
  }
  
  data2plot %>% 
    ggplot(aes(x = variant, y = log_lrt, color = af)) +
    geom_point() +
    labs(
      x = 'Genes',
      y = '-log10(lrt p-value)'
    ) +
    geom_hline(yintercept = c(0,limit)) +
    geom_text_repel(aes(label = labels), box.padding = 1.1,
                    max.overlaps = 100) +
    theme_classic() +
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())
}

manhplot.pval = function(data, limit = 6, unk = TRUE) {
  
  # Draws a manhattan plot from a pyseer output
  # limit marks the threshold by which the funciton will add labels
  # unk: wether you want the Unknowns to be plotted or not
  
  data2plot = data %>% 
    arrange(variant) %>% 
    mutate(log_lrt = -log10(`lrt-pvalue`),
           log_pval = -log10(`filter-pvalue`)) %>% 
    arrange(desc(af)) %>% 
    mutate(variant_single = str_split(variant, pattern = '~~~') %>%  map_chr(., 1)) %>% 
    mutate(variant_2 = case_when( variant_single == 'group_5967' ~ 'Group_5967',
                                  grepl("group_",variant_single) ~ 'Unknown',
                                  TRUE ~ variant_single),
           labels = case_when(log_pval > limit ~ variant_2),
           variant = factor(variant, levels = variant))
  
  if (unk == FALSE){
    data2plot = data2plot %>% 
      mutate(labels = case_when(labels == 'Unknown' ~ '',
                                TRUE ~ labels))
  } else if (unk == TRUE) {
    data2plot = data2plot
  }
  
  data2plot %>% 
    ggplot(aes(x = variant, y = log_pval, color = af)) +
    geom_point() +
    labs(
      x = 'Genes',
      y = '-log10 (p-value)'
    ) +
    geom_hline(yintercept = c(0,limit)) +
    geom_text_repel(aes(label = labels), box.padding = 1.1,
                    max.overlaps = 100) +
    theme_classic() +
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())
}

# plot datasets

### lrt ####
# ALL

manhplot.lrt(all_nobio, limit = 3, unk = F)

dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'manhattan_pyseer_lrt_all_nobio.pdf'),
             height = 10, width = 14, useDingbats = FALSE)


manhplot.lrt(all_worm, limit = 3, unk = F)

dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'manhattan_pyseer_lrt_all_worm.pdf'),
             height = 10, width = 14, useDingbats = FALSE)


# AUS

manhplot.lrt(aus_nobio, limit = 3, unk = F)

dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'manhattan_pyseer_lrt_AUS_nobio.pdf'),
             height = 10, width = 14, useDingbats = FALSE)


manhplot.lrt(aus_worm, limit = 3, unk = F)

dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'manhattan_pyseer_lrt_AUS_worm.pdf'),
             height = 10, width = 14, useDingbats = FALSE)

# ECOREF

manhplot.lrt(ecoref_nobio, limit = 3, unk = F)

dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'manhattan_pyseer_lrt_ECOREF_nobio.pdf'),
             height = 10, width = 14, useDingbats = FALSE)


manhplot.lrt(ecoref_worm, limit = 3, unk = F)

dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'manhattan_pyseer_lrt_ECOREF_worm.pdf'),
             height = 10, width = 14, useDingbats = FALSE)



### pval ####
# ALL

manhplot.pval(all_nobio, limit = 6, unk = F)

dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'manhattan_pyseer_pval_all_nobio.pdf'),
             height = 10, width = 14, useDingbats = FALSE)


manhplot.pval(all_worm, limit = 6, unk = F)

dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'manhattan_pyseer_pval_all_worm.pdf'),
             height = 10, width = 14, useDingbats = FALSE)


# AUS

manhplot.pval(aus_nobio, limit = 6, unk = F)

dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'manhattan_pyseer_pval_AUS_nobio.pdf'),
             height = 10, width = 14, useDingbats = FALSE)


manhplot.pval(aus_worm, limit = 6, unk = F)

dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'manhattan_pyseer_pval_AUS_worm.pdf'),
             height = 10, width = 14, useDingbats = FALSE)

# ECOREF

manhplot.pval(ecoref_nobio, limit = 6, unk = F)

dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'manhattan_pyseer_pval_ECOREF_nobio.pdf'),
             height = 10, width = 14, useDingbats = FALSE)


manhplot.pval(ecoref_worm, limit = 6, unk = F)

dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'manhattan_pyseer_pval_ECOREF_worm.pdf'),
             height = 10, width = 14, useDingbats = FALSE)










