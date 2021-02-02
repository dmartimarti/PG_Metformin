
### libraries ####
library(tidyverse) # master library to deal with data frames
library(readxl) # read xlsx or xls files
library(ggrepel) # ggplot add-on, to plot names that don't collapse in same position
library(FactoMineR) # for PCA
library(factoextra) # for PCA
library(here) # usefull to save plots in folders given a root
library(viridis) # color palette package
library(plotly)
library(broom)
library(gtools)
library(openxlsx)
library(skimr)
library(gridExtra)



theme_set(theme_light())


# after having merged all datasets with the python script named worm_image_merge.py,
# read the files having the data and names, and start calculating stats

# # # # # # # # # # # # # # # 
### Getting data in shape ####
# # # # # # # # # # # # # # # 

# data from imaging scores
worm_data = read_delim("raw_data/worm_imaging_total.csv", "\t",  trim_ws = TRUE)
# rename variables with underscores
names(worm_data) = gsub(" ", "_", names(worm_data))

### DATA FILTERING ###
# probability of single worm > 0.7
# size between 4500 and 12000
data_filt_full = worm_data %>%
  filter(Predicted_Class == 'SingleWorm') %>%
  filter(Probability_of_SingleWorm > 0.7,
         Size_in_pixels > 4500 & Size_in_pixels < 12000) %>% 
  filter(Well != 'ron') %>% 
  mutate(Replicate = as.factor(Replicate)) %>%
  # remove all histograms
  select(-(Histogram_of_Intensity_1:Histogram_of_Intensity_63), -(Histogram_of_Intensity_in_neighborhood_0:Histogram_of_Intensity_in_neighborhood_63))


names(worm_data)



## metadata
strain_db = read_excel("raw_data/strain_db.xlsx")
Natural_isolates_layout = read_excel("raw_data/Natural isolates layout.xlsx", 
                                     sheet = "list_long") 
biofilm_strain_annotation = read_excel("biofilm_strain_annotation.xlsx") %>% select(-`...4`) %>% rename(ID = Strain)


# merge and filter Ev experiments
strain_db = strain_db %>%
  left_join(Natural_isolates_layout) %>% 
  left_join(biofilm_strain_annotation)

# join everything to have names, IDs and all that stuff in one table
data_filt = data_filt %>%
  left_join(strain_db) %>%
  select(PG, Well, Metf, ID, Strainname, B_phenotype, Broadphenotype, everything())


write.csv(data_filt, here('exploration', 'raw_data_filtered.csv'), quote = F, row.names = F)

# save statistics list
list_of_datasets = list('data_filtered' = data_filt)

write.xlsx(list_of_datasets, here('exploration', 'raw_data_filtered.xlsx'), colNames = T, rowNames = F) 



# skim data
skim(worm_data)

# test
data_filt %>%
  filter(Strainname == 'OP50')




### As we have done a 4th replicate, let's save the original filtered list in a new variable that won't be changed
# and I'll filter plates 7 and 8 from everything as they were not included into the 4th replicates (all duplicates)

data_filt_full = data_filt

data_filt = data_filt %>% 
  filter(!(PG %in% c('PG7', 'PG8')))

# # # # # # # # # # # #
### Quality checks ####

# some quality checks
# density plot of the probability of single worm by PG and rep
data_filt %>%
  ggplot(aes(Probability_of_SingleWorm, fill = PG, colour = PG)) +
  geom_density(alpha = 0.1) +
  # facet_wrap(vars(Replicate,PG), nrow = 6) +
  facet_grid(vars(Replicate),vars(PG)) +
  theme_classic()

dev.copy2pdf(device = cairo_pdf,
             file = here('exploration', 'SingleWorm_prob_density.pdf'),
             width = 11, height = 7, useDingbats = FALSE)

# density plot of the probability of single worm by PG and rep
data_filt %>%
  group_by(PG, Replicate) %>%
  # summarise(Mean = mean(Probability_of_SingleWorm, na.rm = TRUE),
  #           SD = sd(Probability_of_SingleWorm, na.rm =  TRUE)) %>%
  ggplot(aes(y = Probability_of_SingleWorm, x = PG, fill = Replicate, group = interaction(Replicate, PG))) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 1) +
  geom_point(position = position_jitterdodge(), alpha = 0.1) +
  theme_classic()


dev.copy2pdf(device = cairo_pdf,
             file = here('exploration', 'SingleWorm_prob_boxplot.pdf'),
             width = 7, height = 7, useDingbats = FALSE)

## having seen this, I would filter every worm that has a prob < 0.7 or so


# well order from 1 to 12, with package naturalsort
lvls = naturalsort::naturalsort(unique(data_filt$Well))

# plot number of worms per plate and replicate
data_filt %>%
  # filter(Probability_of_SingleWorm > 0.65) %>%
  filter(PG == 'PG6') %>%
  mutate(Well = factor(Well, levels = lvls)) %>%
  group_by(Well, PG, Replicate, Metf) %>%
  summarise(N = n()) %>%
  ggplot(aes(x = Replicate, y = N, colour = Metf, fill = Metf)) +
  # geom_point(position = position_jitterdodge(), color = 'black') +
  geom_bar(stat = 'identity', width=.5, position = "dodge") +
  geom_hline(yintercept = 3, linetype = 'dashed', colour = 'grey50') +
  facet_wrap(vars(Well), ncol = 12) +
  theme_classic()


# iterate and save number of worms per plate, well, condition, and replicate
plates = unique(data_filt$PG)

for (plate in plates){
  p = data_filt %>%
    # filter(Probability_of_SingleWorm > 0.65) %>%
    filter(PG == plate) %>%
    mutate(Well = factor(Well, levels = lvls)) %>%
    group_by(Well, PG, Replicate, Metf) %>%
    summarise(N = n()) %>%
    ggplot(aes(x = Replicate, y = N, fill = Metf)) +
    # geom_point(position = position_jitterdodge(), color = 'black') +
    geom_bar(stat = 'identity', width=.85, position = "dodge", color = 'black') +
    geom_hline(yintercept = 3, linetype = 'dashed', colour = 'grey50') +
    scale_fill_manual(values = c('#F5B607', '#0B97D4')) +
    facet_wrap(vars(Well), ncol = 12) +
    theme_classic()
  ggsave(file = here('exploration/Number_of_worms', paste0('N_worms_',plate,'_replicates.pdf')),
         width = 11, height = 7)
}

# number of worms per condition and plate (reps joined)
for (plate in plates){
  p = data_filt %>%
    # filter(Probability_of_SingleWorm > 0.6) %>%
    filter(PG == plate) %>%
    mutate(Well = factor(Well, levels = lvls)) %>%
    group_by(Well, PG, Metf) %>%
    summarise(N = n()) %>%
    ggplot(aes(x = Metf, y = N, fill = Metf)) +
    # geom_point(position = position_jitterdodge()) +
    geom_bar(stat = 'identity', width=.85, position = "dodge", color = 'black') +
    # ylim(0,25) +
    geom_hline(yintercept = 3, linetype = 'dashed', colour = 'grey50') +
    scale_fill_manual(values = c('#F5B607', '#0B97D4')) +
    facet_wrap(vars(Well), ncol = 12) +
    theme_classic()
  ggsave(file = here('exploration/Number_of_worms', paste0('N_worms_',plate,'_condition.pdf')),
         width = 11, height = 9)
}


### EXPLORATION ####

# # # # # # # # # # # #
### Mean Intensity ####
# # # # # # # # # # # #

strains = c('OP50', 'Nissle1917', 'NRG857C', 'MG1655', 'ECOR-02a')
asdf = c('DE-COMM-3445', 'HM-347', 'ECOR-67', 'NILS21')
# strains for Jen's master thesis plot
strains = c('MG1655', 'HM-366', 'DE-COMM-4977', 'OP50', 'NILS19')

# plot Mean intensity by itself (it looks good!)
data_filt %>% 
  # filter(Probability_of_SingleWorm > 0.7) %>%
  filter(Strainname %in% strains) %>%
  # filter(ID %in% IDs) %>% 
  mutate(Strainname = factor(Strainname, levels = strains)) %>%
  ggplot(aes(y = Mean_Intensity, x = Strainname, fill = Metf)) +
  geom_boxplot(outlier.colour = 'red') +
  geom_point(position = position_jitterdodge(jitter.width = 0.3)) +
  scale_fill_manual(values = c('#F5B607', '#0B97D4')) +
  # facet_wrap(~Replicate) +
  theme_light() 

ggsave(file = here('exploration', 'Mean_Intensity_test.pdf'),
       width = 11, height = 7)

# normalise against mean or median
norm = data_filt %>%
  filter(PG == 'PG1') %>%
  # filter(Probability_of_SingleWorm > 0.7) %>%
  group_by(Metf, Replicate) %>%
  summarise(Mean = mean(Mean_Intensity),
            SD = sd(Mean_Intensity),
            Median = median(Mean_Intensity)) %>%
  filter(Metf == '0mM') %>%
  ungroup %>%
  select(!Metf)

# represent the data with different normalisations
# it looks it's a bit better with median, but it's almost the same
data_filt %>% 
  filter(PG == 'PG1') %>%
  # filter(Replicate == 1) %>%
  filter(Strainname %in% strains) %>%
  # filter(Probability_of_SingleWorm > 0.7) %>%
  left_join(norm) %>%
  mutate(median_intensity = Mean_Intensity - Median,
         mean_intensity = Mean_Intensity - Mean) %>%
  ggplot(aes(y = median_intensity, x = Strainname, fill = Metf)) +
  geom_boxplot(outlier.colour = 'red') +
  geom_point(position = position_jitterdodge(jitter.width = 0.3)) +
  # facet_wrap(~Replicate) +
  scale_fill_manual(values = c('#F5B607', '#0B97D4')) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust  = 1))

ggsave(file = here('exploration', 'Mean_Intensity_MedianCorr_test.pdf'),
       width = 11, height = 7)


## Normalise
# plot everything
norm = data_filt %>%
  group_by(PG,Metf, Replicate) %>%
  summarise(Mean = mean(Mean_Intensity),
            SD = sd(Mean_Intensity),
            Median = median(Mean_Intensity)) %>%
  filter(Metf == '0mM') %>%
  ungroup %>%
  select(!Metf)


# loop for median
for (plate in plates){
  data_filt %>% 
    filter(PG == plate) %>%
    # filter(Replicate == 1) %>%
    # filter(Strainname %in% strains) %>%
    left_join(norm) %>%
    mutate(median_intensity = Mean_Intensity - Median,
           mean_intensity = Mean_Intensity - Mean) %>%
    arrange(desc(median_intensity)) %>%
    mutate(Strainname = factor(Strainname)) %>%
    ggplot(aes(y = median_intensity, x = ID, fill = Metf)) +
    geom_boxplot(outlier.colour = 'red') +
    geom_point(position = position_jitterdodge(jitter.width = 0.3), alpha = 0.3) +
    ylim(-1500, 12500) +
    scale_fill_manual(values = c('#F5B607', '#0B97D4')) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90, hjust  = 1))
  
  ggsave(file = here('exploration/Mean_Intensity_Median', paste0('Mean_Intensity_',plate,'_MedianCorrected.pdf')),
         width = 19, height = 7)
}

# loop without correction
for (plate in plates){
  data_filt %>% 
    filter(PG == plate) %>%
    # filter(Replicate == 1) %>%
    # filter(Strainname %in% strains) %>%
    mutate(Strainname = factor(Strainname)) %>%
    ggplot(aes(y = Mean_Intensity, x = ID, fill = Metf)) +
    geom_boxplot(outlier.colour = 'red') +
    geom_point(position = position_jitterdodge(jitter.width = 0.3), alpha = 0.3) +
    ylim(0, 15000) +
    scale_fill_manual(values = c('#F5B607', '#0B97D4')) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90, hjust  = 1))
  
  ggsave(file = here('exploration/Mean_Intensity', paste0('Mean_Intensity_',plate,'.pdf')),
         width = 19, height = 7)
}


## plot every strain in a single page

strains = metadata %>%
  select(ID) %>%
  drop_na(ID) %>%
  t %>% as.character %>% unique

plot_list = list()
for (i in 1:length(strains)){
  p = data_filt %>% 
    filter(ID == strains[i]) %>%
    # filter(Replicate == 1) %>%
    # filter(Probability_of_SingleWorm > 0.85) %>%
    left_join(norm) %>%
    mutate(median_intensity = Mean_Intensity - Median,
           mean_intensity = Mean_Intensity - Mean) %>%
    ggplot(aes(y = Mean_Intensity, x = ID, fill = Metf, group = Metf)) +
    geom_boxplot(outlier.colour = 'pink') +
    geom_point(aes(color = Replicate, size = Size_in_pixels), position = position_jitterdodge(jitter.width = 0.3)) +
    # facet_wrap(~PG) +
    scale_fill_manual(values = c('#F5B607', '#0B97D4')) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 0, hjust  = 1))
  
  plot_list[[i]] = p
}

# save multiple plots into one big file
ggsave(file = here('exploration', 'Mean_Intensity_ALL_STRAINS.pdf'), marrangeGrob(grobs = plot_list, nrow = 1, ncol = 1),
       width = 11, height = 7)




## PLOT WITHOUT Mean_Intensity_in_neighborhood
plot_list = list()
for (i in 1:length(strains)){
  p = data_filt %>% 
    filter(ID == strains[i]) %>%
    # filter(Replicate == 1) %>%
    # filter(Probability_of_SingleWorm > 0.85) %>%
    left_join(norm) %>%
    mutate(median_intensity = Mean_Intensity - Median,
           mean_intensity = Mean_Intensity - Mean) %>%
    ggplot(aes(y = Mean_Intensity - Mean_Intensity_in_neighborhood, x = ID, fill = Metf, group = Metf)) +
    geom_boxplot(outlier.colour = 'pink') +
    geom_point(aes(color = Replicate, size = Size_in_pixels), position = position_jitterdodge(jitter.width = 0.3)) +
    # facet_wrap(~PG) +
    scale_fill_manual(values = c('#F5B607', '#0B97D4')) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 0, hjust  = 1))
  
  plot_list[[i]] = p
}

# save multiple plots into one big file
ggsave(file = here('exploration', 'Mean_Intensity_ALL_STRAINS_noNeighborhoodIntensity.pdf'), marrangeGrob(grobs = plot_list, nrow = 1, ncol = 1),
       width = 11, height = 7)







# data_filt %>% 
#   # select(-Histogram_of_Intensity_0:-Histogram_of_Intensity_63) %>% 
#   select(Skewness_of_Intensity,Minimum_intensity:Bounding_Box_Minimum_0, Mean_Intensity_in_neighborhood, phylogroup) %>%
#   filter(!phylogroup %in% c('albertii', 'Unknown', 'Non Escherichia', 'fergusonii')) %>% 
#   ggpairs(mapping = aes(colour = phylogroup))
# 

data_filt %>% 
  drop_na(phylogroup) %>% 
  filter(!phylogroup %in% c('albertii', 'Unknown', 'Non Escherichia', 'fergusonii',
                            'cladeI', 'cladeII', 'cladeIII', 'cladeIV', 'cladeV')) %>% 
  filter(!(Probability_of_SingleWorm < 0.95)) %>%
  # filter(phylogroup == 'B2') %>% 
  ggplot(aes(x = Size_in_pixels, y = Mean_Intensity, colour = phylogroup, 
             fill = phylogroup, group = phylogroup)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(se = T, method = 'loess', color = 'black', fill = 'grey50') +
  facet_wrap(~Metf+phylogroup, nrow = 2) +
  scale_fill_viridis(discrete = T) + 
  scale_color_viridis(discrete = T)

ggsave(file = here('exploration', 'MeanIntensity_by_size_phylogroup.pdf'),
       width = 15, height = 7)


# # # # # # # # # # # # # #
### duplicated strains ####
# # # # # # # # # # # # # #


dups = metadata %>%
  group_by(ID) %>%
  summarise(N = n()) %>%
  arrange(desc(N)) %>%
  drop_na(ID) %>%
  filter(N > 1) %>% select(ID) %>% t %>% as.character

plot_list = list()
for (i in 1:length(dups)){
  p = data_filt %>% 
    filter(ID == dups[i]) %>%
    # filter(Replicate == 1) %>%
    # filter(Probability_of_SingleWorm > 0.85) %>%
    left_join(norm) %>%
    mutate(median_intensity = Mean_Intensity - Median,
           mean_intensity = Mean_Intensity - Mean) %>%
    ggplot(aes(y = Mean_Intensity, x = ID, fill = Metf, group = Metf)) +
    geom_boxplot(outlier.colour = 'pink') +
    geom_point(aes(color = Replicate, size = Size_in_pixels), position = position_jitterdodge(jitter.width = 0.3)) +
    facet_wrap(~PG) +
    scale_fill_manual(values = c('#F5B607', '#0B97D4')) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 0, hjust  = 1))
  
  plot_list[[i]] = p
}

# save multiple plots into one big file
ggsave(file = here('exploration', 'Mean_Intensity_duplicates.pdf'), marrangeGrob(grobs = plot_list, nrow = 1, ncol = 1),
       width = 11, height = 7)






# # # # # # # # # # # # # # # #  
### Multi-univariate stats ####
# # # # # # # # # # # # # # # # 




# per plate and ID
stat_res = data_filt %>%
  group_by(PG, ID, Strainname, B_phenotype, Broadphenotype, Well) %>%
  nest() %>%
  mutate(models = map(.x = data, .f = lm, formula = 'Mean_Intensity ~ Metf')) %>%
  # mutate(results = map(.f = exploration, .x = models)) %>%	
  mutate(results = map(.f = tidy, .x = models)) %>%
  dplyr::select(PG, Well, Strainname, Broadphenotype, results) %>%
  unnest(cols = c(results)) %>%
  filter(term != '(Intercept)') %>% 
  ungroup %>% # ungroup before passing the FDR correction or it will correct line by line
  mutate(FDR = p.adjust(p.value, method = 'fdr'),
         FDR_stars = gtools::stars.pval(FDR))

stat_res = data_filt %>%
  group_by(PG, Metf, Well) %>%
  summarise(Mean = mean(Mean_Intensity, na.rm = TRUE)) %>%
  pivot_wider(names_from = c(Metf), values_from = Mean, names_prefix = 'Met_') %>%
  mutate(FC = Met_50mM/Met_0mM,
         log2FC = log2(FC)) %>%
  select(PG, Well, FC, log2FC, Met_0mM, Met_50mM) %>%
  left_join(stat_res) %>%
  ungroup %>%
  select(PG, Well, ID, Strainname, Broadphenotype, B_phenotype, FC, log2FC,Met_0mM, Met_50mM, estimate:FDR_stars) %>%
  left_join(data_filt %>% 
              select(PG, ID, Strainname, Broadphenotype, Well, phylogroup) %>% unique)


# no duplicates
stat_res2 = data_filt %>%
  group_by(ID, Strainname, B_phenotype, Broadphenotype) %>%
  nest() %>%
  mutate(models = map(.x = data, .f = lm, formula = 'Mean_Intensity ~ Metf')) %>%
  # mutate(results = map(.f = exploration, .x = models)) %>%	
  mutate(results = map(.f = tidy, .x = models)) %>%
  dplyr::select(ID, Strainname, Broadphenotype, results) %>%
  unnest(cols = c(results)) %>%
  filter(term != '(Intercept)') %>% 
  ungroup %>% # ungroup before passing the FDR correction or it will correct line by line
  mutate(FDR = p.adjust(p.value, method = 'fdr'),
         FDR_stars = gtools::stars.pval(FDR))

stat_res2 = data_filt %>%
  group_by(ID, Strainname, Metf) %>%
  summarise(Mean = mean(Mean_Intensity, na.rm = TRUE)) %>%
  pivot_wider(names_from = c(Metf), values_from = Mean, names_prefix = 'Met_') %>%
  mutate(FC = Met_50mM/Met_0mM,
         log2FC = log2(FC)) %>%
  select(ID, Strainname, FC, log2FC, Met_0mM, Met_50mM) %>%
  left_join(stat_res2) %>%
  ungroup %>%
  select(ID, Strainname, Broadphenotype, B_phenotype, FC, log2FC, Met_0mM, Met_50mM, estimate:FDR_stars) %>%
  left_join(data_filt %>% 
              select(ID, Strainname, Broadphenotype, phylogroup) %>% unique)




# metadata for stats
tables = c('Stats_per_plate', 'Stats_per_strain')
description = c('Statistics with individual duplicates, not corrected', 
                'stats with duplicates merged into single IDs/Strain names')
stats_metadata = data.frame(tables, description)

# save statistics list
list_of_datasets = list('metadata' = stats_metadata, 'Stats_per_plate' = stat_res, 'Stats_per_strain' = stat_res2)

write.xlsx(list_of_datasets, here('exploration', 'worm_imaging_stats.xlsx'), colNames = T, rowNames = F) 


# check some duplicated strains
stat_res %>% 
  filter(ID == dups[1])






# # # # # # # # # # # # # # #  
### phylotypes intensity ####
# # # # # # # # # # # # # # # 

## plot from raw worm intensity

data_filt %>% 
  drop_na(phylogroup) %>%
  filter(!phylogroup %in% c('albertii', 'Unknown', 'Non Escherichia', 'fergusonii')) %>%
  filter(Probability_of_SingleWorm > 0.75) %>%
  left_join(norm) %>%
  ggplot(aes(y = Mean_Intensity, x = phylogroup, fill = phylogroup, group = phylogroup)) +
  geom_boxplot(outlier.colour = 'pink') +
  geom_point(position = position_jitterdodge(jitter.width = 0.4), alpha = 0.1) +
  facet_wrap(~Metf) +
  # scale_fill_manual(values = c('#F5B607', '#0B97D4')) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust  = 1))


ggsave(file = here('exploration', 'Mean_intensity_phylogroup.pdf'),
       width = 15, height = 7)


# statistical tests for phylogroups
test = data_filt %>%
  filter(Metf == '50mM') %>%
  filter(phylogroup %in% c('B1', 'B2', 'A', 'C', 'D', 'E', 'F', 'G'))

fit = aov(Mean_Intensity ~ phylogroup, data = test)
tidy(TukeyHSD(fit)) %>% 
  mutate(pstars = stars.pval(adj.p.value)) %>%View


# how many strains per phylogroup
data_filt %>% 
  drop_na(phylogroup) %>%
  group_by(phylogroup) %>%
  select(ID, phylogroup) %>%
  unique %>%
  summarise(N = n())



# with duplicates
stat_res %>%
  drop_na(phylogroup) %>%
  filter(!phylogroup %in% c('albertii', 'Unknown', 'Non Escherichia', 'fergusonii')) %>%
  # filter(FDR < 0.05) %>%
  ggplot(aes(y = log2FC, x = phylogroup, fill = phylogroup, group = phylogroup)) +
  geom_boxplot(outlier.colour = 'pink') +
  geom_point(position = position_jitterdodge(jitter.width = 0.4), alpha = 0.5) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust  = 1))


ggsave(file = here('exploration', 'log2FC_phylogroup.pdf'),
       width = 13, height = 7)

# per strain, without duplicates
stat_res2 %>%
  drop_na(phylogroup) %>%
  filter(!phylogroup %in% c('albertii', 'Unknown', 'Non Escherichia', 'fergusonii')) %>%
  # filter(FDR < 0.05) %>%
  ggplot(aes(y = log2FC, x = phylogroup, fill = phylogroup, group = phylogroup)) +
  geom_boxplot(outlier.colour = 'pink') +
  geom_point(position = position_jitterdodge(jitter.width = 0.4), alpha = 0.5) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust  = 1))


ggsave(file = here('exploration', 'log2FC_phylogroup_per_strain.pdf'),
       width = 13, height = 7)







# Empty? ------------------------------------------------------------------

data_filt %>%  
  filter(PG ==  'PG7', Well %in% c( 'H2', 'H3')) %>% 
  ggplot(aes(y = Mean_Intensity, x = Metf)) + geom_boxplot() + geom_point(aes(colour = Replicate)) + facet_wrap(~Well)



empty = c('NT12025', 'NT12065', 'NT12066', # PG1
          'NT12105', 'NT12135', 'NT12139', 'NT12175', 'NT12181', # PG2
          'NT12182', 'NT12183', 'NT12184', 'NT12185', 'NT12186', # PG2
          'NT12188', 'NT12190',# PG2
          'NT12548', 'NT12571', 'NT12568' # PG6
)

sim_str = c('NT12102','NT12105','NT12117','NT12117','NT12044',
            'NT12131','NT12375','NT12062','NT12065','NT12175',
            'NT12180','NT12190','NT12228','NT12437','NT12548','NT12571','NT12568')


data_filt %>%  
  filter(ID %in% sim_str) %>% 
  filter(PG %in% c('PG1', 'PG2', 'PG6')) %>% 
  ggplot(aes(y = Mean_Intensity, x = Metf)) + geom_boxplot(aes(fill = Metf)) + 
  geom_point(aes(colour = Replicate)) + 
  facet_wrap(~ID)




# strain filtering --------------------------------------------------------

# calculate the median of the means at 0mM
stat_res %>% 
  # group_by(PG, Well, ID) %>% 
  summarise(Median = median(Met_0mM),
            Mean = mean(Met_0mM),
            SD = sd(Met_0mM))

# cut off = median + 1.5(sd)
cutoff = 2550+(1473 * 1.5)

abnormal = stat_res %>% filter(Met_0mM > cutoff) %>% 
  arrange(Met_0mM)


# save statistics list
list_of_datasets = list('abnormal' = abnormal)

write.xlsx(list_of_datasets, here('exploration', 'abrnormal_strains.xlsx'), colNames = T, rowNames = F) 






# Replicate difference ----------------------------------------------------

replicate_test = data_filt %>% 
  group_by(PG, Well, Metf, ID, Replicate) %>%
  summarise(Mean = mean(Mean_Intensity, na.rm = TRUE)) %>% 
  drop_na(ID) %>% 
  # ungroup %>% 
  pivot_wider(names_from = Replicate, names_prefix = 'rep_', values_from = c(Mean)) %>% 
  mutate(Mean = mean(c(rep_1, rep_2, rep_3), na.rm = T),
         SD = sd(c(rep_1, rep_2, rep_3), na.rm = T)) %>% 
  ungroup %>% 
  arrange(PG, ID) %>% 
  filter(SD > 3000)


# save statistics list
list_of_datasets = list('Replicates' = replicate_test)

write.xlsx(list_of_datasets, here('exploration', 'replicate_test.xlsx'), colNames = T, rowNames = F) 





# Jen plots ---------------------------------------------------------------


data.sum.reps = data_filt %>% 
  group_by(ID, Strainname, Replicate, Metf) %>% 
  summarise(Mean = mean(Mean_Intensity, na.rm = T)) 


data.sum.reps = data.sum.reps %>% 
  filter(Metf == '50mM') %>% 
  pivot_wider(names_from = Replicate, values_from = Mean, names_prefix = 'rep_') %>% 
  drop_na(rep_1, rep_2, rep_3) 



data.sum.reps %>% 
  ggplot(aes(x = rep_2, rep_3)) +
  geom_point() +
  geom_smooth(method = 'lm')

ggsave(file = here('exploration', 'cor_rep2_3.pdf'),
       width = 10, height = 7)




rep1_2 = glance(lm(rep_1 ~ rep_2, data = data.sum.reps))
rep1_3 = glance(lm(rep_1 ~ rep_3, data = data.sum.reps))
rep2_3 = glance(lm(rep_2 ~ rep_3, data = data.sum.reps))

rep_corr = bind_rows(rep1_2, rep1_3, rep2_3) %>% 
  mutate(Contrast = c('rep1_2', 'rep1_3', 'rep2_3'), .before = r.squared)

# save statistics list
list_of_datasets = list('Rep_corr' = rep_corr)

write.xlsx(list_of_datasets, here('exploration', 'replicate_correlation.xlsx'), colNames = T, rowNames = F) 

# 3D plot
library(plotly)
fig <- plot_ly(data = data.sum.reps, x = ~rep_1, y = ~rep_2, z = ~rep_3,
               marker = list(size = 6,
                             color = 'rgba(255, 182, 193, .9)',
                             line = list(color = 'rgba(152, 0, 0, .8)',
                                         width = 2)))

fig


data_filt %>% 
  group_by(ID, Strainname, Metf, phylogroup) %>% 
  summarise(Mean = mean(Mean_Intensity, na.rm = T),
            SD = sd(Mean_Intensity, na.rm = T),
            SEM = SD/sqrt(n())) %>% 
  filter(Metf == '50mM') %>% 
  arrange(desc(Mean)) %>% 
  ungroup %>% 
  mutate(ID = factor(ID, levels = ID)) %>% 
  ggplot(aes(x = ID, y = Mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM)) +
  theme_classic() +
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())

ggsave(file = here('exploration', 'Mean_Int_50mM_SEM.pdf'),
       width = 15, height = 12)

# mean AUC sum from strains, from Jen (with PG6 fixed)
library(readxl)
mean_sumAUC <- read_excel("D:/MRC_Postdoc/Pangenomic/Correlation/200709_mean_sumAUC_PG16.xlsx")


# join both datasets
data_filt %>% 
  group_by(ID, Strainname, Metf, phylogroup) %>% 
  summarise(Mean = mean(Mean_Intensity, na.rm = T),
            SD = sd(Mean_Intensity, na.rm = T),
            SEM = SD/sqrt(n())) %>% 
  filter(Metf == '50mM') %>% 
  arrange(desc(Mean)) %>% 
  ungroup %>% 
  mutate(ID = factor(ID, levels = ID)) %>% 
  rename(Mean_Intensity = Mean,
         SD_Intensity = SD,
         SEM_Intensity = SEM) %>% 
  left_join(mean_sumAUC %>% select(-ID) %>% 
              rename(ID = Strain)) %>% 
  drop_na(Mean, Mean_Intensity) %>% 
  distinct(ID, .keep_all = T) %>% 
  arrange(desc(Mean)) %>% 
  mutate(ID = factor(ID, levels = ID),
         Mean = Mean * 3900) %>% 
  ggplot(aes(x = ID, y = Mean)) +
  geom_point(colour = 'grey') +
  geom_point(aes(x = ID, y = Mean_Intensity), colour = 'red') +
  scale_y_continuous('Worm mean intensity\n (50mM)', sec.axis = sec_axis(~./3000, name = 'Mean resistance')) +
  theme_classic() +
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())




# join both datasets
merge = data_filt %>% 
  group_by(ID, Strainname, Metf, phylogroup) %>% 
  summarise(Mean = mean(Mean_Intensity, na.rm = T),
            SD = sd(Mean_Intensity, na.rm = T),
            SEM = SD/sqrt(n())) %>% 
  filter(Metf == '50mM') %>% 
  arrange(desc(Mean)) %>% 
  rename(Mean_Intensity = Mean,
         SD_Intensity = SD,
         SEM_Intensity = SEM) %>% 
  left_join(mean_sumAUC %>% select(-ID) %>% 
              rename(ID = Strain)) %>% 
  drop_na(Mean, Mean_Intensity, Broadphenotype) %>% 
  distinct(ID, .keep_all = T) %>% 
  filter(Broadphenotype != 'Unknown') 

merge %>% 
  ggplot(aes(x = Mean, y = Mean_Intensity)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  labs(y = 'Worm mean intensity \n (50mM)',
       x = 'Bacterial resistance score') +
  facet_wrap(~Broadphenotype)

ggsave(file = here('exploration', 'Correlations.pdf'), width = 15, height = 7)

corr_stats = merge %>% 
  group_by(Broadphenotype) %>% 
  nest() %>% 
  mutate(fit = map(data, ~lm(Mean_Intensity ~ Mean, data = .x)),
         results = map(fit, glance)) %>% 
  select(-data, -fit) %>% 
  unnest(results) 

# save statistics list
list_of_datasets = list('Phenotype' = corr_stats)

write.xlsx(list_of_datasets, here('exploration', 'phenotype_correlation.xlsx'), colNames = T, rowNames = F) 






# Resistant strains -------------------------------------------------------


res_strains = c('NT12025','NT12051','NT12008','NT12388','NT12045','NT12074','NT12010','NT12353',
                'NT12424','NT12292','NT12321','NT12105','NT12320','NT12360','NT12062','NT12066',
                'NT12447','NT12364','NT12059','NT12538','NT12060','NT12056','NT12294', 'NT12009')

# plot Mean intensity by itself (it looks good!)
data_filt %>% 
  # filter(Probability_of_SingleWorm > 0.7) %>%
  filter(ID %in% res_strains) %>%
  # filter(ID %in% IDs) %>% 
  # mutate(ID = factor(ID, levels = ID)) %>%
  ggplot(aes(y = Mean_Intensity, x = ID, fill = Metf)) +
  geom_boxplot(outlier.colour = 'red') +
  geom_point(position = position_jitterdodge(jitter.width = 0.3)) +
  scale_fill_manual(values = c('#F5B607', '#0B97D4')) +
  # facet_wrap(~Replicate) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file = here('exploration', 'Resistant_strains.pdf'),
       width = 13, height = 7)










# Saving the project ------------------------------------------------------


test.strains = c("NT12038","NT12039","NT12030","NT12031","NT12032","NT12034","NT12035","NT12036",
                 "NT12037","NT12050","NT12051","NT12042","NT12044","NT12045","NT12046","NT12047",
                 "NT12048","NT12049","NT12062","NT12063","NT12054","NT12056","NT12057","NT12058",
                 "NT12059","NT12074","NT12066","NT12067","NT12068","NT12069","NT12070","NT12071",
                 "NT12072","NT12073","NT12086","NT12087","NT12078","NT12079","NT12080","NT12081",
                 "NT12082","NT12084","NT12085","NT12011","NT12002","NT12004","NT12005","NT12006",
                 "NT12007","NT12008","NT12009",'NT12043','NT12055','NT12060','NT12061','NT12075',
                 'NT12083','NT12003','NT12010')



data_filt %>% filter(ID %in% test.strains) %>% 
  arrange(Mean_Intensity) %>% 
  # unique(ID) %>% 
  # mutate(ID = factor(ID, levels = ID)) %>%
  ggplot(aes(y = Mean_Intensity, x = ID, fill = Metf)) +
  geom_boxplot(outlier.colour = 'red') +
  geom_point(position = position_jitterdodge(jitter.width = 0.3)) +
  scale_fill_manual(values = c('#F5B607', '#0B97D4')) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave(file = here('exploration', 'selected_strains.pdf'),
       width = 19, height = 10)



sensitive = c('NT12009','NT12054','NT12068','NT12010','NT12008','NT12050')

resistant = c('NT12032','NT12066','NT12067','NT12062','NT12070')

data_filt %>% filter(ID %in% resistant) %>% 
  arrange(Mean_Intensity) %>% 
  # unique(ID) %>% 
  # mutate(ID = factor(ID, levels = ID)) %>%
  ggplot(aes(y = Mean_Intensity, x = ID, fill = Metf)) +
  geom_boxplot(outlier.colour = 'red') +
  geom_point(position = position_jitterdodge(jitter.width = 0.3)) +
  scale_fill_manual(values = c('#F5B607', '#0B97D4')) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



data.test.sum = data_filt %>% filter(ID %in% test.strains) %>% 
  group_by(ID, Metf) %>% 
  summarise(Mean = mean(Mean_Intensity, na.rm = TRUE),
            SD = sd(Mean_Intensity, na.rm = T)) %>% 
  pivot_wider(names_from = Metf, values_from = c(Mean, SD)) %>% 
  mutate(FC = Mean_50mM / Mean_0mM)


# save statistics list
list_of_datasets = list('worm.test.exploration' = data.test.sum)

write.xlsx(list_of_datasets, here('exploration', 'test_exploration.xlsx'), colNames = T, rowNames = F) 


# Strains > 4000 plots ----------------------------------------------------

res_str =  stat_res2 %>% filter(Met_0mM > 4000) %>% distinct(ID, Met_0mM) %>% select(ID) %>% t %>% as.character()
res_str2 =  stat_res %>% 
  filter(PG %in% c('PG1','PG2','PG3','PG4','PG5','PG6')) %>% 
  filter(Met_0mM > 4000) %>% 
  distinct(ID, Met_0mM) %>% 
  select(ID) %>% t %>% as.character()


## plot every strain in a single page

plot_list = list()
for (i in 1:length(res_str)){
  p = data_filt %>% 
    filter(ID == res_str[i]) %>%
    # filter(Replicate == 1) %>%
    ggplot(aes(y = Mean_Intensity - Mean_Intensity_in_neighborhood, x = ID, fill = Metf, group = Metf)) +
    geom_boxplot(outlier.colour = 'pink') +
    geom_point(aes(color = Replicate, size = Size_in_pixels), position = position_jitterdodge(jitter.width = 0.3)) +
    # facet_wrap(~PG) +
    scale_fill_manual(values = c('#F5B607', '#0B97D4')) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 0, hjust  = 1))
  
  plot_list[[i]] = p
}

# save multiple plots into one big file
ggsave(file = here('exploration', 'Mean_Intensity_resistant_4000_noBackground.pdf'), marrangeGrob(grobs = plot_list, nrow = 1, ncol = 1),
       width = 11, height = 7)




# Biofilm exploration -----------------------------------------------------

stat_res2 %>% 
  mutate(B_phenotype = as.factor(B_phenotype)) %>% 
  filter(B_phenotype != 'contaminate') %>% 
  filter(Met_0mM > 4000) %>%
  ggplot(aes(x = B_phenotype, y = Met_0mM, fill = B_phenotype)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.3))

stat_res2 %>% 
  mutate(B_phenotype = as.factor(B_phenotype)) %>% 
  filter(B_phenotype != 'contaminate') %>% 
  ggplot(aes(y = Met_0mM, x = B_phenotype, fill = B_phenotype)) +
  geom_violin(scale = "width") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.05)




