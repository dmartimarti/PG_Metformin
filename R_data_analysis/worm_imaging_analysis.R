
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
data_filt = worm_data %>%
  filter(Predicted_Class == 'SingleWorm') %>%
  filter(Probability_of_SingleWorm > 0.7,
         Size_in_pixels > 4500 & Size_in_pixels < 12000) %>% 
  mutate(Replicate = as.factor(Replicate)) %>%
  # remove all histograms
  select(-(Histogram_of_Intensity_1:Histogram_of_Intensity_63), -(Histogram_of_Intensity_in_neighborhood_0:Histogram_of_Intensity_in_neighborhood_63))





## metadata
strain_db = read_excel("raw_data/strain_db.xlsx")
Natural_isolates_layout <- read_excel("raw_data/Natural isolates layout.xlsx", 
                                      sheet = "list_long") 

# merge and filter Ev experiments
strain_db = strain_db %>%
  left_join(Natural_isolates_layout) 

# join everything to have names, IDs and all that stuff in one table
data_filt = data_filt %>%
  left_join(strain_db) %>%
  select(PG, Well, Metf, ID, Strainname, Broadphenotype, everything())

# skim data
skim(data_filt)

# test
data_filt %>%
  filter(Strainname == 'OP50')

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
             file = here('summary', 'SingleWorm_prob_density.pdf'),
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
             file = here('summary', 'SingleWorm_prob_boxplot.pdf'),
             width = 7, height = 7, useDingbats = FALSE)

## having seen this, I would filter every worm that has a prob < 0.7 or so


# well order from 1 to 12, with package naturalsort
lvls = naturalsort::naturalsort(unique(data_filt$Well))

# plot number of worms per plate and replicate
data_filt %>%
  # filter(Probability_of_SingleWorm > 0.65) %>%
  filter(PG == 'PG7') %>%
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
  ggsave(file = here('summary/Number_of_worms', paste0('N_worms_',plate,'_replicates.pdf')),
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
  ggsave(file = here('summary/Number_of_worms', paste0('N_worms_',plate,'_condition.pdf')),
         width = 11, height = 9)
}


### EXPLORATION ####

# # # # # # # # # # # #
### Mean Intensity ####
# # # # # # # # # # # #

strains = c('OP50', 'Nissle1917', 'NRG857C', 'MG1655', 'ECOR-02a')
asdf = c('DE-COMM-3445', 'HM-347', 'ECOR-67', 'NILS21')

# plot Mean intensity by itself (it looks good!)
data_filt %>% 
  # filter(Probability_of_SingleWorm > 0.7) %>%
  filter(Strainname %in% strains) %>%
  mutate(Strainname = factor(Strainname, levels = strains)) %>%
  ggplot(aes(y = Mean_Intensity, x = Strainname, fill = Metf)) +
  geom_boxplot(outlier.colour = 'red') +
  geom_point(position = position_jitterdodge(jitter.width = 0.3)) +
  scale_fill_manual(values = c('#F5B607', '#0B97D4')) +
  # facet_wrap(~Replicate) +
  theme_light() 

ggsave(file = here('summary', 'Mean_Intensity_test.pdf'),
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

ggsave(file = here('summary', 'Mean_Intensity_MedianCorr_test.pdf'),
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
  
  ggsave(file = here('summary/Mean_Intensity_Median', paste0('Mean_Intensity_',plate,'_MedianCorrected.pdf')),
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
  
  ggsave(file = here('summary/Mean_Intensity', paste0('Mean_Intensity_',plate,'.pdf')),
         width = 19, height = 7)
}


## plot every strain in a single page

strains = Natural_isolates_layout %>%
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
ggsave(file = here('summary', 'Mean_Intensity_ALL_STRAINS.pdf'), marrangeGrob(grobs = plot_list, nrow = 1, ncol = 1),
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

ggsave(file = here('summary', 'MeanIntensity_by_size_phylogroup.pdf'),
       width = 15, height = 7)


# # # # # # # # # # # # # #
### duplicated strains ####
# # # # # # # # # # # # # #


dups = Natural_isolates_layout %>%
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
    filter(Probability_of_SingleWorm > 0.85) %>%
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
ggsave(file = here('summary', 'Mean_Intensity_duplicates_prob085.pdf'), marrangeGrob(grobs = plot_list, nrow = 1, ncol = 1),
       width = 11, height = 7)






# # # # # # # # # # # # # # # #  
### Multi-univariate stats ####
# # # # # # # # # # # # # # # # 




# per plate and ID
stat_res = data_filt %>%
  group_by(PG, ID, Strainname, Broadphenotype, Well) %>%
  nest() %>%
  mutate(models = map(.x = data, .f = lm, formula = 'Mean_Intensity ~ Metf')) %>%
  mutate(results = map(.f = summary, .x = models)) %>%	
  mutate(results = map(.f = tidy, .x = results)) %>%
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
  select(PG, Well, ID, Strainname, Broadphenotype, FC, log2FC,Met_0mM, Met_50mM, estimate:FDR_stars) %>%
  left_join(data_filt %>% 
              select(PG, ID, Strainname, Broadphenotype, Well, phylogroup) %>% unique)
  

# no duplicates
stat_res2 = data_filt %>%
  group_by(ID, Strainname, Broadphenotype) %>%
  nest() %>%
  mutate(models = map(.x = data, .f = lm, formula = 'Mean_Intensity ~ Metf')) %>%
  mutate(results = map(.f = summary, .x = models)) %>%	
  mutate(results = map(.f = tidy, .x = results)) %>%
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
  select(ID, Strainname, Broadphenotype, FC, log2FC, Met_0mM, Met_50mM, estimate:FDR_stars) %>%
  left_join(data_filt %>% 
              select(ID, Strainname, Broadphenotype, phylogroup) %>% unique)




# metadata for stats
tables = c('Stats_per_plate', 'Stats_per_strain')
description = c('Statistics with individual duplicates, not corrected', 
                'stats with duplicates merged into single IDs/Strain names')
stats_metadata = data.frame(tables, description)

# save statistics list
list_of_datasets = list('metadata' = stats_metadata, 'Stats_per_plate' = stat_res, 'Stats_per_strain' = stat_res2)

write.xlsx(list_of_datasets, here('Summary', 'worm_imaging_stats.xlsx'), colNames = T, rowNames = F) 


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


ggsave(file = here('summary', 'Mean_intensity_phylogroup.pdf'),
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
  

ggsave(file = here('summary', 'log2FC_phylogroup.pdf'),
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


ggsave(file = here('summary', 'log2FC_phylogroup_per_strain.pdf'),
       width = 13, height = 7)









