# This script analyses the results from pyseer and plots them as 
# Manhattan plots along the population abundance per gene


# libraries ---------------------------------------------------------------



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
library(glue)

theme_set(theme_cowplot(14))


# metadata = read_excel("D:/MRC_Postdoc/Pangenomic/pangenome_analysis/metadata/MAIN_metadata.xlsx", 
#                       sheet = "metadata")

metadata = read_excel("~/Documents/MRC_postdoc/Pangenomic/metadata/MAIN_metadata.xlsx", 
                      sheet = "metadata")


# phylogroups -------------------------------------------------------------


# summarise phylogroups 
phylogr = metadata %>%
  filter(Discard == 'No') %>% 
  filter(phylogroup != 'cladeI') %>% 
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
  filter(phylo_corrected != 'cladeI') %>% 
  mutate(Total = sum(N),
         Fraction = round((N/Total)*100,1)) %>%
  group_by(phylo_corrected) %>% 
  summarise(Fraction = sum(Fraction)) %>% 
  mutate(y_label_pos = Fraction + 3) %>% 
  ggplot(aes(x = phylo_corrected, y = Fraction, fill = phylo_corrected)) +
  geom_bar(stat = 'identity', color = 'black') +
  geom_text(aes(y = y_label_pos, label = Fraction), vjust = 1.6, size = 3.5) +
  labs(x = 'Phylogroups', y = '% of total') +
  theme_cowplot(18) +
  guides(fill = 'none') +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1))


dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'phylogrous_dist_barplot.pdf'),
             height = 6, width = 8, useDingbats = FALSE)





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
  drop_na(phylo_corrected, Broadphenotype) %>% 
  group_by(Broadphenotype) %>% 
  mutate(Total = sum(N),
         Fraction = round((N/Total)*100,1),
         y_label_pos = Fraction + 6) %>%
  filter(Broadphenotype != 'Unknown') %>% 
  ggplot(aes(x = phylo_corrected, y = Fraction, fill = phylo_corrected)) +
  geom_bar(stat = 'identity', color = 'black') +
  geom_text(aes(y = y_label_pos, label = Fraction), vjust = 1.6, size = 3.5) +
  labs(x = 'phylogroups', y = '% of total') +
  facet_wrap(~Broadphenotype, ncol = 3) +
  theme_cowplot(18) +
  guides(fill='none') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'phylogrous_phenotype_dist_barplot.pdf'),
             height = 6, width = 10, useDingbats = FALSE)




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


# biofilm producers
biofilm_prod = read_delim("pyseer_output/results/worm_phenotype_normal2superbio.tsv", 
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


# super bio producers

manhplot.lrt(biofilm_prod, limit = 3, unk = F)

dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'manhattan_pyseer_lrt_normal2superbio.pdf'),
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


# super bio producers

manhplot.pval(biofilm_prod, limit = 6, unk = F)

dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'manhattan_pyseer_pval_normal2superbio.pdf'),
             height = 10, width = 14, useDingbats = FALSE)




### extract genes for enrichment with String DB

# helper function to extract genes from pyseer
get_genes = function(pyseer_out = all_nobio, 
                     mode = 'lrt', # can be either pval or lrt
                     thres = 3 # -log10 threshold to filter
                     ){
  
  if (mode == 'lrt'){
    temp = pyseer_out %>% 
      mutate(log_lrt = -log10(`lrt-pvalue`),
             log_pval = -log10(`filter-pvalue`)) %>% 
      filter(log_lrt > thres) %>% 
      filter(!str_detect(variant, 'group')) %>% 
      select(variant) %>% 
      mutate(variant = str_split(variant, '~~~')) %>% 
      unnest(variant) %>% 
      mutate(variant = case_when(str_detect(variant, '_') ~ str_sub(variant, end = -3),
                                 TRUE ~ variant),
             variant = str_replace(variant, '\\d','')) %>% 
      distinct(variant)
  } else if (mode == 'pval') {
    temp = pyseer_out %>% 
      mutate(log_lrt = -log10(`lrt-pvalue`),
             log_pval = -log10(`filter-pvalue`)) %>% 
      filter(log_pval > thres) %>% 
      filter(!str_detect(variant, 'group')) %>% 
      select(variant) %>% 
      mutate(variant = str_split(variant, '~~~')) %>% 
      unnest(variant) %>% 
      mutate(variant = case_when(str_detect(variant, '_') ~ str_sub(variant, end = -3),
                                 TRUE ~ variant),
             variant = str_replace(variant, '\\d','')) %>% 
      distinct(variant)
  } else {
    cat("You need to introduce lrt or pval as a mode! \n")
  }
  
  return(temp)

}


#### save datasets ####

pval_thres = 6
lrt_thres = 3

# ALL
get_genes(all_nobio, mode = 'lrt', thres = lrt_thres) %>% 
  write_delim( here('pyseer_output/gene_lists',glue('all_nobio_lrt_{lrt_thres}.txt')), 
               col_names = F)
get_genes(all_nobio, mode = 'pval', thres = pval_thres) %>% 
  write_delim( here('pyseer_output/gene_lists',glue('all_nobio_pval_{pval_thres}.txt')), 
               col_names = F)


get_genes(all_worm, mode = 'lrt', thres = lrt_thres) %>% 
  write_delim( here('pyseer_output/gene_lists',glue('all_worm_lrt_{lrt_thres}.txt')), 
               col_names = F)
get_genes(all_worm, mode = 'pval', thres = pval_thres) %>% 
  write_delim( here('pyseer_output/gene_lists',glue('all_worm_pval_{pval_thres}.txt')), 
               col_names = F)

# AUS
get_genes(aus_nobio, mode = 'lrt', thres = lrt_thres) %>% 
  write_delim( here('pyseer_output/gene_lists',glue('aus_nobio_lrt_{lrt_thres}.txt')), 
               col_names = F)
get_genes(aus_nobio, mode = 'pval', thres = pval_thres) %>% 
  write_delim( here('pyseer_output/gene_lists',glue('aus_nobio_pval_{pval_thres}.txt')), 
               col_names = F)


get_genes(aus_worm, mode = 'lrt', thres = lrt_thres) %>% 
  write_delim( here('pyseer_output/gene_lists',glue('aus_worm_lrt_{lrt_thres}.txt')), 
               col_names = F)
get_genes(aus_worm, mode = 'pval', thres = pval_thres) %>% 
  write_delim( here('pyseer_output/gene_lists',glue('aus_worm_pval_{pval_thres}.txt')), 
               col_names = F)


# pyseer bact ------------------------------------------------------------------


### read data ####

# ALL no bio
all_nobio_bact = read_delim("pyseer_output/results/bact/bact_phenotype_no_biofilm.tsv", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)

# ALL
all_bact = read_delim("pyseer_output/results/bact/bact_phenotype_ALL.tsv", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)


# superbio producers
supebio_bact = read_delim("pyseer_output/results/bact/bact_phenotype_normal2superbio.tsv", 
                      "\t", escape_double = FALSE, trim_ws = TRUE)



#### lrt plot ####


# super bio producers

manhplot.lrt(supebio_bact, limit = 3, unk = F)

dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'manhattan_pyseer_BACT_lrt_normal2superbio.pdf'),
             height = 10, width = 14, useDingbats = FALSE)



##### pval plot ####

# super bio producers

manhplot.pval(supebio_bact, limit = 6, unk = F)

dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'manhattan_pyseer_BACT_pval_normal2superbio.pdf'),
             height = 10, width = 14, useDingbats = FALSE)











# Gene families -----------------------------------------------------------



# load the gene presence matrix

gene_pa =read_csv("gene_presence_absence.csv")

# names(gene_pa)

gene_pa = gene_pa %>% 
  select(-Annotation, -`Non-unique Gene name`, gene = Gene, everything())

gene_bin = gene_pa %>% 
  select(-gene) %>% 
  mutate_all(~replace(.,!is.na(.), 1)) %>% 
  mutate_all(~replace(., is.na(.), 0)) %>% 
  # mutate_all(~replace(., is.character(.), 1)) %>% 
  mutate_all(~as.integer(.)) 

genes = gene_pa$gene


rw_sums = rowSums(gene_bin)


gene_bin = gene_bin %>% 
  mutate(total = rw_sums, .before = `100`) %>% 
  mutate(gene = genes, .before = total)



gene_bin %>% 
  ggplot(aes(x = total)) +
  geom_histogram(position = 'identity',
                 bins=70,
                 fill = 'grey60', 
                 color = 'black') +
  labs(y = 'Gene count',
       x = 'Number of genomes with a specific gene') +
  theme_cowplot(20)

dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'gene_counts_per_genome.pdf'),
             height = 5, width = 6, useDingbats = FALSE)








# how many genes are included in each gene family

genes_split = str_split(genes, '~~~')

gene_fam_lengths = c()
for(element in genes_split){
  # print(length(element))
  gene_fam_lengths = c(gene_fam_lengths ,length(element))
}

gene_fam_lengths


gene_fam_df = tibble(gene = genes, lengths = gene_fam_lengths)

gene_fam_df %>% 
  arrange(desc(lengths))


gene_fam_df %>% 
  filter(lengths > 1) %>% 
  ggplot(aes(x = lengths)) +
  geom_histogram(
    stat = 'count'
    ) +
  # scale_x_discrete(breaks=c(seq(2,17,1))) +
  scale_x_continuous(limits=c(1, 18),
                     breaks=c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)) +
  labs(
    x = 'Number of genes in family',
    y = 'Count'
  )

dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'gene_per_genefamily.pdf'),
             height = 8, width = 10, useDingbats = FALSE)





gene_fam_df %>% 
  # filter(lengths > 1) %>% 
  ggplot(aes(x = lengths)) +
  geom_histogram(
    stat = 'count'
  ) +
  # scale_x_discrete(breaks=c(seq(2,17,1))) +
  # scale_x_continuous(limits=c(0, 18),
  #                    breaks=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)) +
  labs(
    x = 'Number of genes in family',
    y = 'Count'
  ) + 
  ggforce::facet_zoom(x = lengths > 2, 
                      xlim = c(2.2,17),
                      ylim = c(0,2300))

dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'gene_per_genefamily_zoom.pdf'),
             height = 8, width = 13, useDingbats = FALSE)




### what are the genes that belong to the shell 


gene_bin %>% 
  select(gene, total) %>% 
  mutate(per = total / 750) %>% 
  filter(per < 0.95 & per >= 0.15) %>% 
  write_csv('shell_genes.csv')




## is there any correlation between genes within families and gene presence? 

gene_fam_df = gene_fam_df %>% 
  left_join(gene_bin %>% select(gene, total)) %>% 
  mutate(per = total / 750) %>% 
  mutate(class = case_when( per >= 0.99 ~ 'core',
                            per < 0.99 & per >= 0.95 ~ 'soft_core',
                            per < 0.95 & per >= 0.15 ~ 'shell',
                            per < 0.15 ~ 'cloud'))




gene_fam_df %>%
  ggplot(aes(x = per, y = lengths)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  ggpubr::stat_cor(
    p.accuracy = 0.001, r.accuracy = 0.01
  ) +
  labs(
    x = 'Gene presence in pangenome',
    y = 'Genes within a gene family'
  ) +
  theme_cowplot(14)


dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'genePA_length_corr.pdf'),
             height = 8, width = 13, useDingbats = FALSE)


gene_fam_df %>%
  filter(class != 'core') %>%
  ggplot(aes(x = per, y = lengths, color = class)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = 'lm', color = 'black') +
  ggpubr::stat_cor(
    p.accuracy = 0.001, r.accuracy = 0.01) +
  labs(
    x = 'Gene presence in pangenome',
    y = 'Genes within a gene family'
  ) +
  theme_cowplot(14) + 
  facet_wrap(~class, scales = 'free_x')


dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'genePA_length_corr_class.pdf'),
             height = 8, width = 13, useDingbats = FALSE)





gene_fam_df %>% 
  mutate(known = case_when(str_detect(gene, 'group') == TRUE ~ 'Unknown',
                           TRUE ~ 'Known')) %>% 
  group_by(known) %>% 
  count() %>% 
  ggplot(aes(x=known, y = n)) +
  geom_col(aes(fill = known), color = 'black') +
  guides(fill = 'none') + 
  scale_fill_manual(values = c('#2869EB', '#EBD420')) +
  labs(
    x = 'Gene source',
    y = 'Number of gene families'
  ) + 
  theme_cowplot(19)


dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'known_unknown_genes.pdf'),
             height = 7, width = 6, useDingbats = FALSE)




gene_fam_df %>% 
  mutate(known = case_when(str_detect(gene, 'group') == TRUE ~ 'Unknown',
                           TRUE ~ 'Known'),
         class = as_factor(class),
         known = as_factor(known)) %>% 
  # group_by(known, class) %>% 
  # count() %>% 
  ggplot(aes(x = known, fill = class)) +
  geom_histogram(stat = 'count', 
                 position = 'dodge',
                 color = 'black') + 
  labs(
    x = 'Gene source',
    y = 'Number of gene families'
  ) + 
  theme_cowplot(17)


dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'known_unknown_genes_perClass.pdf'),
             height = 8, width = 9, useDingbats = FALSE)



# proportion of genes known or unknown by class
gene_fam_df %>% 
  mutate(known = case_when(str_detect(gene, 'group') == TRUE ~ 'Unknown',
                           TRUE ~ 'Known'),
         class = as_factor(class),
         known = as_factor(known)) %>% 
  group_by(class, known) %>% 
  count() %>% 
  group_by(class) %>% 
  mutate(tot = sum(n)) %>% 
  ungroup %>% 
  mutate(prop = (n / tot) * 100) %>% 
  ggplot(aes(x = known, y = prop, fill = class)) +
  geom_col(color = 'black',
           position = 'dodge') +
  geom_text(aes(label = round(prop,2)), 
            position = position_dodge(width = 0.9),
            vjust=-0.25, 
            size = 4) + 
  ylim(0,100) +
  # guides(fill = 'none') + 
  labs(
    x = 'Gene source',
    y = '% of gene families'
  ) + 
  theme_cowplot(19)

dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'known_unknown_genes_perClass_prop.pdf'),
             height = 7, width = 7, useDingbats = FALSE)







gene_fam_df %>% 
  filter(str_detect(gene, 'bcs'))







# PCA of gene P/A ---------------------------------------------------------

gene_bin2 = gene_bin %>% 
  select(-total)

# transpose the matrix to have gene names as columns
gene_bin2 = gene_bin2 %>% 
  pivot_longer(cols = where(is.numeric), 
               values_to = 'presence', names_to = 'genome') %>% 
  pivot_wider(names_from = gene, values_from = presence)


pca_fit <- gene_bin2 %>% 
  select(where(is.numeric)) %>% # retain only numeric columns
  prcomp(scale = F)




pca_plot = pca_fit %>%
  broom::augment(gene_bin2) %>% # add original dataset back in
  select(genome, .fittedPC1:.fittedPC5) %>% 
  left_join(metadata %>% 
              mutate(genome = str_sub(fasta, start = 1, end = -7), .before = ID) %>% 
              select(genome, Broadphenotype, phylogroup))



pca_plot %>% 
  drop_na(phylogroup) %>% 
  filter(phylogroup != 'cladeI') %>% 
  filter(phylogroup != 'E or cladeI') %>% 
  ggplot(aes(.fittedPC1, .fittedPC2, color = phylogroup, fill = phylogroup)) + 
  geom_point(size = 2, alpha = 0.7) +
  # stat_ellipse() +
  stat_ellipse(level=0.95, geom = 'polygon', alpha = 0.6) +
  theme_half_open(12) 

dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'PCA_gene_PA.pdf'),
             height = 8, width = 9, useDingbats = FALSE)






