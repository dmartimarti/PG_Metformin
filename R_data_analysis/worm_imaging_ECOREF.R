
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


# names(worm_data)



## metadata
strain_db = read_excel("raw_data/strain_db.xlsx")
Natural_isolates_layout = read_excel("raw_data/Natural isolates layout.xlsx", 
                                     sheet = "list_long") %>% 
  filter(PG %in%  c('PG1','PG2','PG3','PG4','PG5','PG6'))

biofilm_strain_annotation = read_excel("20909_biofilm_strain_annotation.xlsx", sheet = 'rep4') %>% rename(ID = Strain) %>% 
  select(-`...7`,-`...8`,-`...9`) %>% 
  rename(biofilm_0mM = `0mM_annotation`,
         biofilm_50mM = `50mM_annotation`)


# merge and filter Ev experiments
metadata = strain_db %>%
  left_join(Natural_isolates_layout) %>% 
  filter(!(is.na(PG) == TRUE)) %>% 
  left_join(biofilm_strain_annotation)


# join everything to have names, IDs and all that stuff in one table
data_filt = data_filt %>%
  left_join(metadata) %>%
  select(PG, Well, Metf, ID, Strainname, B_phenotype, biofilm_0mM, biofilm_50mM, Broadphenotype, everything())


list_of_datasets = list('MAIN_metadata' = metadata)
write.xlsx(list_of_datasets, here('exploration', 'MAIN_metadata.xlsx'), colNames = T, rowNames = F) 

# test
data_filt %>%
  filter(Strainname == 'OP50')



### As we have done a 4th replicate, let's save the original filtered list in a new variable that won't be changed
# and I'll filter plates 7 and 8 from everything as they were not included into the 4th replicates (all duplicates)

data_filt_full = data_filt

data_filt = data_filt %>% 
  filter(!(PG %in% c('PG7', 'PG8')))




# # # # # # # # # # # # # #  
### Outlier detection  ####
# # # # # # # # # # # # # # 


# bunch of functions to detect outliers

# z-score
isnt_out_z <- function(x, thres = 3, na.rm = TRUE) {
  abs(x - mean(x, na.rm = na.rm)) <= thres * sd(x, na.rm = na.rm)
}

# z-score with MAD
isnt_out_mad <- function(x, thres = 3, na.rm = TRUE) {
  abs(x - median(x, na.rm = na.rm)) <= thres * mad(x, na.rm = na.rm)
}

# Tukey fences
isnt_out_tukey <- function(x, k = 1.5, na.rm = TRUE) {
  quar <- quantile(x, probs = c(0.25, 0.75), na.rm = na.rm)
  iqr <- diff(quar)
  
  (quar[1] - k * iqr <= x) & (x <= quar[2] + k * iqr)
}

# Mahalanobis detection (MULTIVARIATE)
maha_dist <- . %>% select_if(is.numeric) %>%
  mahalanobis(center = colMeans(.), cov = cov(.))

isnt_out_maha <- function(tbl, isnt_out_f, ...) {
  tbl %>% maha_dist() %>% isnt_out_f(...)
}


isnt_out_funs <- funs(
  z = isnt_out_z,
  mad = isnt_out_mad,
  tukey = isnt_out_tukey
)

# calculate outliers by 3 different methods
data_filt_outs = data_filt %>% 
  group_by(PG,Well,Metf,ID,Strainname,B_phenotype,biofilm_0mM,biofilm_50mM,Broadphenotype) %>% 
  mutate(z_out = isnt_out_z(Mean_Intensity),
         mad_out = isnt_out_mad(Mean_Intensity),
         tukey_out = isnt_out_tukey(Mean_Intensity)) %>% 
  ungroup

# check how many datapoints we have
sum(data_filt_outs$tukey_out)

data_filt_mad = data_filt_outs %>% filter(mad_out == TRUE)




# WELL IT SEEMS IT WORKS BETTER!
# save the dataset as data_filt and continue with next analyses


# save data in csv and xlsx format
write.csv(data_filt_mad, here('exploration', 'raw_data_filtered.csv'), quote = F, row.names = F)

# save statistics list
list_of_datasets = list('data_filtered' = data_filt_mad)
write.xlsx(list_of_datasets, here('exploration', 'raw_data_filtered.xlsx'), colNames = T, rowNames = F) 



# some quality checks
# density plot of the probability of single worm by PG and rep
data_filt_mad %>%
  ggplot(aes(Probability_of_SingleWorm, fill = PG, colour = PG)) +
  geom_density(alpha = 0.1) +
  # facet_wrap(vars(Replicate,PG), nrow = 6) +
  facet_grid(vars(Replicate),vars(PG)) +
  theme_classic()

dev.copy2pdf(device = cairo_pdf,
             file = here('exploration', 'SingleWorm_prob_density.pdf'),
             width = 11, height = 7, useDingbats = FALSE)

# density plot of the probability of single worm by PG and rep
data_filt_mad %>%
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
data_filt_mad %>%
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
plates = unique(data_filt_mad$PG)

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
# asdf = c('DE-COMM-3445', 'HM-347', 'ECOR-67', 'NILS21')
# # strains for Jen's master thesis plot
# strains = c('MG1655', 'HM-366', 'DE-COMM-4977', 'OP50', 'NILS19')
IDs = c('NT12085')
# plot Mean intensity by itself (it looks good!)
data_filt_mad %>% 
  # filter(Probability_of_SingleWorm > 0.7) %>%
  # filter(Strainname %in% strains) %>%
  filter(ID %in% IDs) %>% 
  mutate(Strainname = factor(Strainname, levels = strains)) %>%
  ggplot(aes(y = Mean_Intensity - Mean_Intensity_in_neighborhood , x = Strainname, fill = Metf)) +
  geom_boxplot(outlier.colour = 'red') +
  geom_point(position = position_jitterdodge(jitter.width = 0.3)) +
  scale_fill_manual(values = c('#F5B607', '#0B97D4')) +
  # facet_wrap(~Replicate) +
  theme_light() 

# ggsave(file = here('exploration', 'Mean_Intensity_test.pdf'),
#        width = 11, height = 7)



## plot every strain in a single page

strains = metadata %>%
  select(ID) %>%
  drop_na(ID) %>%
  t %>% as.character %>% unique

plot_list = list()
for (i in 1:length(strains)){
  p = data_filt_mad %>% 
    filter(ID == strains[i]) %>%
    # filter(Replicate == 1) %>%
    # filter(Probability_of_SingleWorm > 0.85) %>%
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
ggsave(file = here('analysis', 'Mean_Intensity_ALL_STRAINS.pdf'), marrangeGrob(grobs = plot_list, nrow = 1, ncol = 1),
       width = 11, height = 7)




## PLOT WITHOUT Mean_Intensity_in_neighborhood
plot_list = list()
for (i in 1:length(strains)){
  p = data_filt_mad %>% 
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



data_filt_mad %>% 
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
### Summary statistics ####
# # # # # # # # # # # # # #


# this piece of code calculates FC per replicate, and then calculate the mean 
# of the FCs, its SD and also the median
FC_means = data_filt_mad %>% 
  group_by(Metf, Replicate, PG, Well, ID, Strainname, biofilm_0mM, biofilm_50mM) %>% 
  summarise(Mean = mean(Mean_Intensity)) %>% 
  pivot_wider(names_from = Metf, values_from = Mean, names_prefix = 'Metf_') %>% 
  mutate(FC = Metf_50mM/Metf_0mM) %>%
  group_by(PG, Well, ID, Strainname, biofilm_0mM, biofilm_50mM)  %>% 
  summarise(Mean_FC = mean(FC, na.rm = TRUE),
            SD_FC = sd(FC, na.rm = TRUE),
            Median_FC = median(FC, na.rm = TRUE))

FC_means %>% 
  arrange(FC_means) %>% 
  ggplot(aes(x = ID, y = Mean_FC)) +
  geom_pointrange(aes(ymin = Mean_FC - SD_FC, ymax = Mean_FC + SD_FC))

# save a list of unique IDs 
FC_means_unique = FC_means %>% 
  ungroup %>% 
  distinct(ID, .keep_all = T)

write_csv(FC_means, here('analysis', 'FC_means.csv'))

write_csv(FC_means_unique, here('analysis', 'FC_means_unique.csv'))



sum.stats = data_filt_mad %>% 
  group_by(Metf, ID, Strainname, biofilm_0mM, biofilm_50mM) %>% 
  summarise(Mean = mean(Mean_Intensity),
            SD = sd(Mean_Intensity)) 


list_of_datasets = list(
  Summ_stats = sum.stats,
  FC_means = FC_means
)

library(openxlsx)

write.xlsx(list_of_datasets, here('analysis', 'Summary_Stats_FC.xlsx'))


# # # # # # # # # # # # # # # #  
### Multi-univariate stats ####
# # # # # # # # # # # # # # # # 


# per plate and ID
stat_res = data_filt_mad %>%
  group_by(PG, ID, Strainname, B_phenotype, biofilm_0mM, biofilm_50mM, Broadphenotype, Well) %>%
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

stat_res = data_filt_mad %>%
  group_by(PG, Metf, Well) %>%
  summarise(Mean = mean(Mean_Intensity, na.rm = TRUE)) %>%
  pivot_wider(names_from = c(Metf), values_from = Mean, names_prefix = 'Met_') %>%
  mutate(FC = Met_50mM/Met_0mM,
         log2FC = log2(FC)) %>%
  select(PG, Well, FC, log2FC, Met_0mM, Met_50mM) %>%
  left_join(stat_res) %>%
  ungroup %>%
  select(PG, Well, ID, Strainname, Broadphenotype, B_phenotype, biofilm_0mM, biofilm_50mM, FC, log2FC,Met_0mM, Met_50mM, estimate:FDR_stars) %>%
  left_join(data_filt_mad %>% 
              select(PG, ID, Strainname, Broadphenotype, Well, phylogroup) %>% unique)


# no duplicates
stat_res2 = data_filt_mad %>%
  group_by(ID, Strainname, B_phenotype,biofilm_0mM, biofilm_50mM,  Broadphenotype) %>%
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

stat_res2 = data_filt_mad %>%
  group_by(ID, Strainname, Metf) %>%
  summarise(Mean = mean(Mean_Intensity, na.rm = TRUE)) %>%
  pivot_wider(names_from = c(Metf), values_from = Mean, names_prefix = 'Met_') %>%
  mutate(FC = Met_50mM/Met_0mM,
         log2FC = log2(FC)) %>%
  select(ID, Strainname, FC, log2FC, Met_0mM, Met_50mM) %>%
  left_join(stat_res2) %>%
  ungroup %>%
  select(ID, Strainname, Broadphenotype, B_phenotype,biofilm_0mM, biofilm_50mM,  FC, log2FC, Met_0mM, Met_50mM, estimate:FDR_stars) %>%
  left_join(data_filt_mad %>% 
              select(ID, Strainname, Broadphenotype, phylogroup) %>% unique)

# metadata for stats
tables = c('Stats_per_plate', 'Stats_per_strain')
description = c('Statistics with individual duplicates, not corrected', 
                'stats with duplicates merged into single IDs/Strain names')
stats_metadata = data.frame(tables, description)

# save statistics list
list_of_datasets = list('metadata' = stats_metadata, 'Stats_per_plate' = stat_res, 'Stats_per_strain' = stat_res2)

write.xlsx(list_of_datasets, here('analysis', 'worm_imaging_stats.xlsx'), colNames = T, rowNames = F) 
write.xlsx(list_of_datasets, here('exploration', 'worm_imaging_stats.xlsx'), colNames = T, rowNames = F) 


# check some duplicated strains
stat_res %>% 
  filter(ID == dups[1])



stat_res2 %>% 
  filter(biofilm_50mM == 'normal', 
         Met_0mM < 4000) %>% 
  drop_na(phylogroup) %>%
  filter(!(phylogroup %in% c('cladeI','Non Escherichia'))) %>% 
  write_csv(here('analysis','worm_imaging_SUBSET.csv'))




# sublibrary --------------------------------------------------------------



# data from imaging scores
sub_library = read_delim("subscreen/worm_imaging_subscreen.csv", "\t",  trim_ws = TRUE)
# rename variables with underscores
names(sub_library) = gsub(" ", "_", names(sub_library))

### DATA FILTERING ###
# probability of single worm > 0.7
# size between 4500 and 12000
sub_library = sub_library %>%
  filter(Predicted_Class == 'SingleWorm') %>%
  filter(Probability_of_SingleWorm > 0.7,
         Size_in_pixels > 4500 & Size_in_pixels < 12000) %>% 
  filter(Well != 'ron') %>% 
  select(-PG, -X1)


# names(worm_data)



## metadata

sub_metadata = read_excel('subscreen/201020_biofilm_annotation_sub_library.xlsx') 

# join everything to have names, IDs and all that stuff in one table
sub_library = sub_library %>%
  left_join(sub_metadata %>% filter(Replicate_bio == 1)) %>%
  left_join(metadata) %>% 
  select(PG, Well, Metf, Strainname, Strain, biofilm_0mM, biofilm_50mM, Broadphenotype, everything()) %>% 
  mutate(Replicate = as.factor(Replicate))

# test
sub_library %>%
  filter(Strainname == 'OP50')



## EXPLORATION PLOTS FOR SUBLIBRARY ##

# total n per condition
sub_library %>% 
  count(Well,ID, Metf) %>% 
  ggplot(aes(x = Metf, y = n, fill = Metf)) +
  geom_histogram(stat= 'identity', colour = 'black') +
  geom_hline(yintercept = 3) +
  facet_wrap(~Well, ncol = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file = here('exploration', 'sublib_numbers_of_worms.pdf'),
       width = 10, height = 9)


# total n per replicate
sub_library %>% 
  count(Well,ID, Replicate, Metf) %>% 
  ggplot(aes(x = Metf, y = n, fill = Replicate)) +
  geom_histogram(stat= 'identity', colour = 'black', position=position_dodge()) +
  facet_wrap(~Well, ncol = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave(file = here('exploration', 'sublib_numbers_of_worms_byreplicate.pdf'),
       width = 10, height = 9)


## Mean brightness per worm

sub_library %>% 
  ggplot(aes(x = Metf, y = Mean_Intensity, fill = Metf)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(), alpha = .4) +
  facet_wrap(~Well, ncol = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust =  1))

ggsave(file = here('exploration', 'sublib_Mean_Intensity.pdf'),
       width = 16, height = 15)


##

sub_strains = unique(sub_library$Strain)

ori_data = data_filt_mad %>% filter(ID %in% sub_strains) %>% 
  select(ID, Metf, Replicate, Mean_Intensity) %>% 
  mutate(Origin = 'original')


sub_data = sub_library %>% select(ID = Strain, Metf, Replicate, Mean_Intensity) %>% 
  mutate(Origin = 'sublibrary')


unique(ori_data$ID)
unique(sub_data$ID)

setdiff(unique(ori_data$ID),unique(sub_data$ID))



sublib_test = ori_data %>% bind_rows(sub_data)

sublib_test %>% unite('condition',Metf,Origin,remove = FALSE) %>% 
  ggplot(aes(x = condition, y = Mean_Intensity, fill = condition)) +
  geom_boxplot() +
  geom_point(alpha = 0.5, aes(colour = Replicate)) +
  geom_vline(xintercept = 2.5) +
  facet_wrap(~ID, ncol = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust =  1))


ggsave(file = here('exploration', 'sublib_Mean_Intensity_COMPARISON.pdf'),
       width = 18, height = 15)



# # # # # # #
## explore biofilm formation in these strains

# distances from sublibrary
sub_dist = sub_metadata %>% 
  # filter(Strain %in% sub_strains) %>% 
  filter(Replicate_bio == 1) %>% 
  rename(biofilm_0mM = Annotation_0mM,
         biofilm_50mM = Annotation_50mM) %>% 
  mutate(biofilm_0mM = case_when(biofilm_0mM == 'normal' | biofilm_0mM == 'NB' ~ 0,
                                 biofilm_0mM == 'biofilm' ~ 1,
                                 biofilm_0mM == 'super_bio' ~ 2),
         biofilm_50mM = case_when(biofilm_50mM == 'normal' ~ 0,
                                  biofilm_50mM == 'biofilm' ~ 1,
                                  biofilm_50mM == 'super_bio' ~ 2),
         distance = biofilm_50mM - biofilm_0mM,
         data = 'sublibrary') %>% 
  arrange(desc(distance)) %>% 
  select(Strain, distance, data)

# distances from original data
ori_dist = metadata %>% 
  filter(ID %in% sub_strains) %>% 
  mutate(biofilm_0mM = case_when(biofilm_0mM == 'normal' | biofilm_0mM == 'NB' ~ 0,
                                 biofilm_0mM == 'biofilm' ~ 1,
                                 biofilm_0mM == 'super_bio' ~ 2),
         biofilm_50mM = case_when(biofilm_50mM == 'normal' ~ 0,
                                  biofilm_50mM == 'biofilm' ~ 1,
                                  biofilm_50mM == 'super_bio' ~ 2),
         distance = biofilm_50mM - biofilm_0mM,
         data = 'original') %>% 
  drop_na(ID) %>% 
  arrange(desc(distance)) %>% 
  select(Strain = ID, distance, data)


dists = rbind(sub_dist,ori_dist)

# plot distances

dists %>% 
  ggplot(aes(x = distance, y = Strain, color = data)) +
  geom_point(aes(shape = data), size = 3, alpha = 0.6) +
  scale_color_manual(values=c("blue", "red"))

ggsave(file = here('exploration', 'sublib_biofilm_dist_COMPARISON.pdf'),
       width = 8, height = 15)



























