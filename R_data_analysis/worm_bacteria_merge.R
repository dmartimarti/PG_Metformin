
# libraries ####

library(tidyverse)
library(openxlsx)
library(here)
library(readxl)
library(colorspace)
library(broom)
library(tidymodels)
library(ggpubr)

theme_set(theme_light())



# Read data ---------------------------------------------------------------

### BACTERIA ###
# unweighted scores
# let's use the second method here, more robust
un.resist = read_excel("D:/MRC_Postdoc/Pangenomic/biolog/PG_growth_biospa/Growth_resistance_summaryStats_NODUPS.xlsx", 
                              sheet = "unweighted") %>% select(-...1, -ID) %>% rename(ID = Strain)

# WORM
# raw data
data_filt = read_csv("D:/MRC_Postdoc/Pangenomic/Worm_imaging/exploration/raw_data_filtered.csv")
# summary stats per well
stat_res = read_excel("D:/MRC_Postdoc/Pangenomic/Worm_imaging/analysis/worm_imaging_stats.xlsx", 
                      sheet = "Stats_per_plate")
# summary stats per strain
stat_res2 = read_excel("D:/MRC_Postdoc/Pangenomic/Worm_imaging/analysis/worm_imaging_stats.xlsx", 
                       sheet = "Stats_per_strain")

# metadata with antiSMASH info
metadata = read_excel("D:/MRC_Postdoc/Pangenomic/Worm_imaging/MAIN_metadata.xlsx", 
                            sheet = "metadata_antismash") %>% 
  select(-...1)

# # # # # # # # #
# correlations ------------------------------------------------------------
# # # # # # # # #


###
# FILTER PHYLOGROUPS THAT ARE NOT E COLI
###
corr_df = stat_res2 %>% 
  left_join(un.resist) %>%
  filter(Met_0mM < 4000) %>%
  filter(Mean < 20) %>% 
  filter(!(phylogroup %in% c('cladeI','cladeII','cladeIII','cladeIV','cladeV','fergusonii', 'albertii', 'Non Escherichia')))




# biofilm 50 mM

corr_df %>%
  # filter(ID != 'NT12177')%>%
  # drop_na(phylogroup) %>% 
  ggplot(aes(x = log2(Mean),y = log2FC)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(y = 'log2 Fold change in worm brightness \n (50 mM / 0 mM)',
       x = 'log2 Bacterial growth score') +
  facet_wrap(~biofilm_50mM)


ggsave(file = here('summary', 'log2FC_AUC_correlation_by_biofilm(50mM).pdf'), width = 130, height = 80, units = 'mm', scale = 2, device = 'pdf')


# biofilm 0 mM

corr_df %>%
  # filter(ID != 'NT12177')%>%
  # drop_na(phylogroup) %>% 
  ggplot(aes(x = log2(Mean),y = log2FC)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(y = 'log2 Fold change in worm brightness \n (50 mM / 0 mM)',
       x = 'log2 Bacterial growth score') +
  facet_wrap(~biofilm_0mM)


ggsave(file = here('summary', 'log2FC_AUC_correlation_by_biofilm(0mM).pdf'), width = 130, height = 80, units = 'mm', scale = 2, device = 'pdf')



# Broadphenotype

corr_df %>%
  # filter(ID != 'NT12177')%>%
  # drop_na(phylogroup) %>% 
  ggplot(aes(x = log2(Mean),y = log2FC)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(y = 'log2 Fold change in worm brightness \n (50 mM / 0 mM)',
       x = 'log2 Bacterial growth score') +
  facet_wrap(~Broadphenotype)


ggsave(file = here('summary', 'log2FC_AUC_correlation_by_Phenotype.pdf'), width = 130, height = 80, units = 'mm', scale = 2, device = 'pdf')


###
# by product
###

products = unique(metadata$product)
products = products[!is.na(products)]

for (prod in products){
  corr_df %>% 
    left_join(metadata) %>% 
    filter(product == prod) %>% 
    distinct(ID, .keep_all = T) %>% 
    filter(Met_0mM < 4000) %>%
    filter(Mean < 20) %>%
    # filter(ID != 'NT12177')%>%
    # drop_na(phylogroup) %>% 
    ggplot(aes(x = log2(Mean),y = log2FC)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    labs(y = 'log2 Fold change in worm brightness \n (50 mM / 0 mM)',
         x = 'log2 Bacterial growth score')
  
  ggsave(file = here('summary', paste0('products/log2FC_AUC_correlation_Products_',prod,'.pdf')), 
         width = 80, height = 80, units = 'mm', scale = 2, device = 'pdf')
}

# facet wrap by biofilm 50 mM
for (prod in products){
  corr_df %>%
    left_join(metadata) %>% 
    filter(product == prod) %>% 
    distinct(ID, .keep_all = T) %>% 
    filter(Met_0mM < 4000) %>%
    filter(Mean < 20) %>%
    # filter(ID != 'NT12177')%>%
    # drop_na(phylogroup) %>% 
    ggplot(aes(x = log2(Mean),y = log2FC)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    labs(y = 'log2 Fold change in worm brightness \n (50 mM / 0 mM)',
         x = 'log2 Bacterial growth score') +
    facet_wrap(~biofilm_50mM)
  
  ggsave(file = here('summary', paste0('products_biofilm50mM/log2FC_AUC_correlation_Products_',prod,'.pdf')), 
         width = 130, height = 80, units = 'mm', scale = 2, device = 'pdf')
}





# univariate stats --------------------------------------------------------

# calculate stats with the different comparisons

# biofilm 50mM
bio50.stats = corr_df %>% 
  mutate(log2Mean = log2(Mean)) %>% 
  group_by(biofilm_50mM) %>% 
  nest() %>% 
  mutate(models = map(.x = data, .f = lm, formula = 'log2FC ~ log2Mean'),
         results = map(.f = glance, .x = models)) %>% 
  select(variable = biofilm_50mM, results) %>% 
  mutate(grouping = 'biofilm_50mM') %>% 
  unnest(cols = c(results)) %>%
  mutate(correlation = sqrt(r.squared)) %>% 
  ungroup %>% # ungroup before passing the FDR correction or it will correct line by line
  mutate(p_stars = gtools::stars.pval(p.value))

# biofilm 0mM
bio0.stats = corr_df %>% 
  mutate(log2Mean = log2(Mean)) %>% 
  group_by(biofilm_0mM) %>% 
  nest() %>% 
  mutate(models = map(.x = data, .f = lm, formula = 'log2FC ~ log2Mean'),
         results = map(.f = glance, .x = models)) %>% 
  select(variable = biofilm_0mM, results) %>% 
  mutate(grouping = 'biofilm_0mM') %>%
  unnest(cols = c(results)) %>%
  mutate(correlation = sqrt(r.squared)) %>% 
  ungroup %>% # ungroup before passing the FDR correction or it will correct line by line
  mutate(p_stars = gtools::stars.pval(p.value))

# phenotype
pheno = corr_df %>% 
  mutate(log2Mean = log2(Mean)) %>% 
  group_by(Broadphenotype) %>% 
  nest() %>% 
  mutate(models = map(.x = data, .f = lm, formula = 'log2FC ~ log2Mean'),
         results = map(.f = glance, .x = models)) %>% 
  select(variable = Broadphenotype, results) %>% 
  mutate(grouping = 'broad phenotype') %>%
  unnest(cols = c(results)) %>%
  mutate(correlation = sqrt(r.squared)) %>%
  ungroup %>% # ungroup before passing the FDR correction or it will correct line by line
  mutate(p_stars = gtools::stars.pval(p.value))

complete.stats = bio50.stats %>% bind_rows(bio0.stats, pheno) %>% select(grouping, variable, everything())

# per product 
prod.df = tibble()
for (prod in products){
  
  temp = corr_df %>% 
    left_join(metadata) %>% 
    filter(product == prod) %>% 
    distinct(ID, .keep_all = T) %>% 
    mutate(log2Mean = log2(Mean)) %>% 
    nest(data = everything()) %>% 
    mutate(models = map(.x = data, .f = lm, formula = 'log2FC ~ log2Mean'),
           results = map(.f = glance, .x = models)) %>% 
    select(results) %>% 
    mutate(grouping = 'products', product = prod) %>%
    unnest(cols = c(results)) %>%
    mutate(correlation = sqrt(r.squared)) %>%
    ungroup %>% # ungroup before passing the FDR correction or it will correct line by line
    mutate(p_stars = gtools::stars.pval(p.value))
  
  prod.df = prod.df %>% bind_rows(temp)

}


# per product and biofilm (50 mM)
prod.bio50.df = tibble()
for (prod in products){
  
  temp = corr_df %>% 
    left_join(metadata) %>% 
    filter(product == prod) %>% 
    distinct(ID, .keep_all = T) %>% 
    mutate(log2Mean = log2(Mean)) %>% 
    group_by(biofilm_50mM) %>% 
    nest() %>% 
    mutate(models = map(.x = data, .f = lm, formula = 'log2FC ~ log2Mean'),
           results = map(.f = glance, .x = models)) %>% 
    select(biofilm_50mM, results) %>% 
    mutate(grouping = 'products', product = prod) %>%
    unnest(cols = c(results)) %>%
    mutate(correlation = sqrt(r.squared)) %>%
    ungroup %>% # ungroup before passing the FDR correction or it will correct line by line
    mutate(p_stars = gtools::stars.pval(p.value))
  
  prod.bio50.df = prod.bio50.df %>% bind_rows(temp)
  
}


# per product and biofilm (0 mM)
prod.bio0.df = tibble()
for (prod in products){
  
  temp = corr_df %>% 
    left_join(metadata) %>% 
    filter(product == prod) %>% 
    distinct(ID, .keep_all = T) %>% 
    mutate(log2Mean = log2(Mean)) %>% 
    group_by(biofilm_0mM) %>% 
    nest() %>% 
    mutate(models = map(.x = data, .f = lm, formula = 'log2FC ~ log2Mean'),
           results = map(.f = glance, .x = models)) %>% 
    select(biofilm_0mM, results) %>% 
    mutate(grouping = 'products', product = prod) %>%
    unnest(cols = c(results)) %>%
    mutate(correlation = sqrt(r.squared)) %>%
    ungroup %>% # ungroup before passing the FDR correction or it will correct line by line
    mutate(p_stars = gtools::stars.pval(p.value))
  
  prod.bio0.df = prod.bio0.df %>% bind_rows(temp)
  
}


# save results in PDF
my_list = list('complete_stats' = complete.stats,
               'antiSMASH_products' = prod.df,
               'antiSMASH_products_biofilm_50mM' = prod.bio50.df,
               'antiSMASH_products_biofilm_0mM' = prod.bio0.df)

write.xlsx(my_list, here('summary', 'correlation_tests.xlsx'), colNames = T)




# # # # # # # # 
# exploration -------------------------------------------------------------
# # # # # # # # 


## FC by biofilm formation
corr_df %>% 
  ggplot(aes(x = biofilm_0mM, y = FC, group = biofilm_0mM, fill = biofilm_0mM)) +
  geom_boxplot(outlier.color = 'red', outlier.shape = 3, outlier.size = 0.5) +
  geom_point(position = position_jitterdodge())

ggsave(file = here('summary', 'worm_FC_by_biofilm0.pdf'), 
       width = 90, height = 80, units = 'mm', scale = 2, device = 'pdf')

corr_df %>% 
  ggplot(aes(x = biofilm_50mM, y = FC, group = biofilm_50mM, fill = biofilm_50mM)) +
  geom_boxplot(outlier.color = 'red', outlier.shape = 3, outlier.size = 0.5) +
  geom_point(position = position_jitterdodge())

ggsave(file = here('summary', 'worm_FC_by_biofilm50.pdf'), 
       width = 90, height = 80, units = 'mm', scale = 2, device = 'pdf')


# by product

prod.stats = corr_df %>% 
  left_join(metadata)

prod.stats %>% 
  distinct(ID, product, .keep_all = T) %>% 
  ggplot(aes(x = product, y = FC, group = product, fill = product)) +
  geom_boxplot(outlier.color = 'red', outlier.shape = 3, outlier.size = 0.5) +
  geom_point(position = position_jitterdodge(), alpha = 0.5) +
  facet_wrap(~biofilm_50mM) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file = here('summary', 'worm_FC_by_product_biofilm50.pdf'), 
       width = 120, height = 80, units = 'mm', scale = 2, device = 'pdf')


prod.stats %>% 
  distinct(ID, product, .keep_all = T) %>% 
  ggplot(aes(x = biofilm_50mM, y = FC, group = biofilm_50mM, fill = biofilm_50mM)) +
  geom_boxplot(outlier.color = 'red', outlier.shape = 3, outlier.size = 0.5) +
  geom_point(position = position_jitterdodge(), alpha = 0.5) +
  facet_wrap(~product) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file = here('summary', 'worm_FC_by_biofilm50_product.pdf'), 
       width = 120, height = 90, units = 'mm', scale = 2, device = 'pdf')




## dataset for detailed phenotype changes
# THIS DATASET INCLUDES ALL WORMS, INCLUDING THE ONES WITH HIGH BRIGHTNESS AT 0 mM
stat_bio = stat_res2 %>% 
  mutate(change = case_when(biofilm_0mM == 'normal' & biofilm_50mM == 'normal' ~ 'no change (normal)',
                            biofilm_0mM == 'NB' & biofilm_50mM == 'normal' ~ 'no change (normal)',
                            biofilm_0mM == 'normal' & biofilm_50mM == 'biofilm' ~ 'biofilm formation',
                            biofilm_0mM == 'normal' & biofilm_50mM == 'super_bio' ~ 'super_bio formation',
                            biofilm_0mM == 'biofilm' & biofilm_50mM == 'biofilm' ~ 'no change (biofilm)',
                            biofilm_0mM == 'biofilm' & biofilm_50mM == 'normal' ~ 'normal from biofilm',
                            biofilm_0mM == 'biofilm' &  biofilm_50mM == 'super_bio' ~ 'super_bio formation (biofilm)',
                            biofilm_0mM == 'super_bio'& biofilm_50mM == 'super_bio' ~ 'no change (super_bio)',
                            biofilm_0mM == 'super_bio'& biofilm_50mM == 'biofilm' ~ 'biofilm from super_bio',
                            biofilm_0mM == 'super_bio'& biofilm_50mM == 'normal' ~ 'normal from super_bio')) %>% 
  drop_na(change)



stat_bio %>% 
  distinct(ID, change) %>% 
  count(change) %>% 
  drop_na(change) %>% 
  ggplot(aes(x = change, y = n, fill = change)) +
  geom_bar(stat = 'identity', colour = 'black') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = 'Biofilm phenotype change (50 mM vs 0 mM)',
       y = 'Number of events in the dataset')

ggsave(file = here('summary', 'biofilm_change.pdf'), 
       width = 110, height = 90, units = 'mm', scale = 2, device = 'pdf')


stat_bio %>% 
  ggplot(aes(x = change, y = Met_50mM, group = change, fill = change)) +
  geom_boxplot(outlier.color = 'red', outlier.shape = 3, outlier.size = 0.5) +
  geom_point(position = position_jitterdodge()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = 'Biofilm phenotype change (50 mM vs 0 mM)',
       y = 'Worm brightness at 50 mM')

ggsave(file = here('summary', 'worm_50mM_(biofilm change).pdf'), 
       width = 90, height = 80, units = 'mm', scale = 2, device = 'pdf')


stat_bio %>% 
  ggplot(aes(x = change, y = Met_0mM, group = change, fill = change)) +
  geom_boxplot(outlier.color = 'red', outlier.shape = 3, outlier.size = 0.5) +
  geom_point(position = position_jitterdodge()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = 'Biofilm phenotype change (50 mM vs 0 mM)',
       y = 'Worm brightness at 0 mM')

ggsave(file = here('summary', 'worm_0mM_(biofilm change).pdf'), 
       width = 90, height = 80, units = 'mm', scale = 2, device = 'pdf')


stat_bio %>% 
  ggplot(aes(x = change, y = FC, group = change, fill = change)) +
  geom_boxplot(outlier.color = 'red', outlier.shape = 3, outlier.size = 0.5) +
  geom_point(position = position_jitterdodge()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = 'Biofilm phenotype change (50 mM vs 0 mM)',
       y = 'Worm brightness at 0 mM')

ggsave(file = here('summary', 'worm_FC_(biofilm change).pdf'), 
       width = 90, height = 80, units = 'mm', scale = 2, device = 'pdf')

# calculate relative distances of phenotype (50mM - 0 mM)
# normal = 0
# biofilm = 1
# super_bio = 2


stat_dist = corr_df %>% 
  mutate(biofilm_0mM = case_when(biofilm_0mM == 'normal' | biofilm_0mM == 'NB' ~ 0,
                                 biofilm_0mM == 'biofilm' ~ 1,
                                 biofilm_0mM == 'super_bio' ~ 2),
         biofilm_50mM = case_when(biofilm_50mM == 'normal' ~ 0,
                                  biofilm_50mM == 'biofilm' ~ 1,
                                  biofilm_50mM == 'super_bio' ~ 2),
         distance = biofilm_50mM - biofilm_0mM) %>% 
  drop_na(ID) %>% 
  arrange(desc(distance))


# comb plot 1 -------------------------------------------------------------

## create the dataset to be ploted with distances
comb = stat_res2 %>% 
  filter(!(phylogroup %in% c('cladeI','cladeII','cladeIII','cladeIV','cladeV','fergusonii', 'albertii', 'Non Escherichia'))) %>% 
  mutate(biofilm_0mM = case_when(biofilm_0mM == 'normal' | biofilm_0mM == 'NB' ~ 0,
                                 biofilm_0mM == 'biofilm' ~ 1,
                                 biofilm_0mM == 'super_bio' ~ 2),
         biofilm_50mM = case_when(biofilm_50mM == 'normal' ~ 0,
                                  biofilm_50mM == 'biofilm' ~ 1,
                                  biofilm_50mM == 'super_bio' ~ 2),
         distance = biofilm_50mM - biofilm_0mM) %>% 
  drop_na(ID) %>% 
  arrange(desc(distance))


 
dplot = comb %>% 
  arrange(desc(distance)) %>% 
  distinct(ID, .keep_all = T) %>% 
  mutate(ID = factor(ID, levels = ID)) %>% 
  ggplot(aes(y = ID, x = distance)) +
  # geom_segment(aes(xend = -2, yend = distance)) +
  geom_point(size = 0.5) +
  theme_classic() + 
  theme(axis.text.y=element_blank(),
        axis.ticks.y = element_blank()) + 
  labs(y = 'Strain',
       x = 'Phenotype distance')


fc.plot = comb %>% 
  arrange(desc(distance)) %>% 
  distinct(ID, .keep_all = T) %>% 
  mutate(ID = factor(ID, levels = ID),
         distance = as.factor(distance)) %>% 
  ggplot(aes(y = ID, x = FC, colour = distance)) +
  # geom_segment(aes(xend = -2, yend = distance)) +
  geom_point() +
  theme_classic() + 
  # xlim(1500, 13000) +
  theme(axis.text.y=element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") +
  # labs(x = 'Worm brightness at 0 mM')
  labs(x = 'Worm brightness FC')


fc.boxplot = comb %>% 
  arrange(desc(distance)) %>% 
  distinct(ID, .keep_all = T) %>% 
  mutate(ID = factor(ID, levels = ID),
         distance = as.factor(distance)) %>% 
  ggplot(aes(y = ID, x = FC, colour = distance, group = distance, fill = distance)) +
  # geom_segment(aes(xend = -2, yend = distance)) +
  geom_boxplot(color = 'black') +
  # geom_point(color = 'black', alpha = .5, position = position_jitterdodge(jitter.width = 0.8)) +
  theme_classic() + 
  # xlim(1500, 13000) +
  theme(axis.text.y=element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
  # labs(x = 'Worm brightness at 0 mM')
  labs(x = 'Worm brightness FC')

ggarrange(dplot, fc.plot, fc.boxplot,
          ncol = 3, nrow = 1, widths = c(0.7,2,2), heights = c(0.9,0.9,0.9))

ggsave(file = here('summary', 'worm_FC_phenotype_distance.pdf'), 
       width = 150, height = 100, units = 'mm', scale = 2, device = 'pdf')





# comb plot 2 -------------------------------------------------------------






dplot = stat_bio %>% 
  left_join(stat_dist %>% select(ID, distance)) %>% 
  arrange(desc(change, distance)) %>%
  distinct(ID, .keep_all = T) %>% 
  mutate(ID = factor(ID, levels = ID)) %>%
  ggplot(aes(y = ID, x = distance)) +
  # geom_segment(aes(xend = -2, yend = distance)) +
  geom_point(size = 0.5) +
  theme_classic() + 
  theme(axis.text.y=element_blank(),
        axis.ticks.y = element_blank()) + 
  labs(y = 'Strain',
       x = 'Phenotype distance')


fc.plot = stat_bio %>% 
  left_join(stat_dist %>% select(ID, distance)) %>% 
  arrange(desc(change, distance)) %>%
  distinct(ID, .keep_all = T) %>% 
  mutate(ID = factor(ID, levels = ID)) %>% 
  ggplot(aes(y = ID, x = FC, colour = change)) +
  # geom_segment(aes(xend = -2, yend = distance)) +
  geom_point() +
  # xlim(1500, 13000) +
  theme_classic() + 
  theme(axis.text.y=element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") +
  # labs(x = 'Worm brightness at 50 mM')
  labs('Worm brightness FC')


fc.boxplot = stat_bio %>% 
  left_join(stat_dist %>% select(ID, distance)) %>%
  arrange(desc(change)) %>%
  distinct(ID, .keep_all = T) %>%
  mutate(ID = factor(ID, levels = ID)) %>%
  ggplot(aes(y = ID, x = FC, group = change, fill = change)) +
  # geom_segment(aes(xend = -2, yend = distance)) +
  geom_boxplot(color = 'black') +
  geom_point(color = 'black', alpha = .5, position = position_jitterdodge(jitter.width = 0.6)) +
  theme_classic() + 
  # xlim(1500, 13000) +
  theme(axis.text.y=element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
  # labs(x = 'Worm brightness at 50 mM')
  labs('Worm brightness FC')

ggarrange(dplot, fc.plot, fc.boxplot,
          ncol = 3, nrow = 1, widths = c(0.7,1.5,2.5), heights = c(0.9,0.9,0.9))

ggsave(file = here('summary', 'worm_FC_phenotype_distance2.pdf'), 
       width = 150, height = 100, units = 'mm', scale = 2, device = 'pdf')






# tSNE --------------------------------------------------------------------
library(Rtsne)

tsne_df = stat_dist %>% 
  left_join(un.resist) %>% 
  distinct(ID, .keep_all = T) %>% 
  select(ID, biofilm_0mM, biofilm_50mM, FC, Met_0mM, Met_50mM, phylogroup, distance, AUC50, AUC100, AUC200, Mean_Bact_metf_0, Mean_Bact_metf_50, Mean_Bact_metf_100, Mean_Bact_metf_200) %>% 
  drop_na(phylogroup) %>% 
  mutate(phylogroup = case_when(phylogroup == 'A' ~ 1,
                                phylogroup == 'B1' ~ 2,
                                phylogroup == 'B2' ~ 3,
                                phylogroup == 'C' ~ 4,
                                phylogroup == 'D' ~ 5,
                                phylogroup == 'E' ~ 6,
                                phylogroup == 'E or cladeI' ~ 7,
                                phylogroup == 'G' ~ 8,
                                phylogroup == 'F' ~ 9,
  )) 

library(ggrepel)
tsne_df

# after playing with the tSNE algorithm, perplexity 10 and theta 3 seem to be optimal to separate some groups in the strains
# I have now a state where I can distinguish some groups, re-run the thing at your own convenience, but paramters 
# down the script will need to be changed accordingly
tsne = Rtsne(tsne_df[2:15], dims = 2, perplexity=10, theta = 0.3, verbose=TRUE, max_iter = 10000, num_threads = 8)

tsne.res = tsne$Y

tsne.res = tsne.res %>% as_tibble() %>% mutate(ID = tsne_df$ID) %>% left_join(stat_dist)


tsne.res 
tsne.res %>% 
  mutate(biofilm_0mM = as.factor(biofilm_0mM),
         biofilm_50mM = as.factor(biofilm_50mM)) %>% 
  ggplot(aes(V1, V2, color = biofilm_50mM)) +
  geom_point(size = 3) 


tsne.res %>% filter(V2 > 40) %>% View

tsne.res %>% filter(V2 > 25) %>% View



ggsave(file = here('summary', 'tSNE.perplex10.biofilm50mM.pdf'), 
       width = 120, height = 100, units = 'mm', scale = 2, device = 'pdf')
















