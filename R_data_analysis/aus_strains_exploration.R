
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


# # # # # # # # # # # # # # # 
### Getting data in shape ####
# # # # # # # # # # # # # # # 

# data from imaging scores
worm_data = read_delim("worm_imaging_AUS.csv", "\t",  trim_ws = TRUE) %>% select(-X1)
# rename variables with underscores
names(worm_data) = gsub(" ", "_", names(worm_data))

### DATA FILTERING ###
# probability of single worm > 0.7
# size between 4500 and 12000
data_filt_full = worm_data %>%
  filter(Predicted_Class == 'SingleWorm') %>%
  filter(Probability_of_SingleWorm > 0.6,Size_in_pixels > 4500 & Size_in_pixels < 12000) %>%
  mutate(Replicate = as.factor(Replicate)) %>% 
  unite('ID', Strain, PG, remove=FALSE) %>% 
  select(-Strain)



## metadata
# strain_db = read_excel("D:/MRC_Postdoc/Pangenomic/Worm_imaging/raw_data/strain_db.xlsx")

# reading only rep1 for now, but I'll need to check if there are differences
biofilm_strain_annotation = read_excel("201201_biofilm_annotation_australian_strains.xlsx", sheet = 'rep1') %>% 
  rename(PG = Plate) %>% 
  mutate(Strain = case_when(Strain == 'OP50' ~ '0',
                            TRUE ~ Strain),
         Strain = as.numeric(Strain)) %>% 
  select(-Replicate) %>% 
  unite('ID', Position, PG, remove=FALSE)


# read the other replicates
# rep 1
annot_1 = read_excel("201201_biofilm_annotation_australian_strains.xlsx", sheet = 'rep1') %>% 
  rename(PG = Plate) %>% 
  mutate(Strain = case_when(Strain == 'OP50' ~ '0',
                            TRUE ~ Strain),
         Strain = as.numeric(Strain)) %>% 
  unite('ID', Position, PG, remove=FALSE)
# rep2
annot_2 = read_excel("201201_biofilm_annotation_australian_strains.xlsx", sheet = 'rep2') %>% 
  rename(PG = Plate) %>% 
  mutate(Strain = case_when(Strain == 'OP50' ~ '0',
                            TRUE ~ Strain),
         Strain = as.numeric(Strain)) %>% 
  unite('ID', Position, PG, remove=FALSE)
# rep3 
annot_3 = read_excel("201201_biofilm_annotation_australian_strains.xlsx", sheet = 'rep3') %>% 
  rename(PG = Plate) %>% 
  mutate(Strain = case_when(Strain == 'OP50' ~ '0',
                            TRUE ~ Strain),
         Strain = as.numeric(Strain)) %>% 
  unite('ID', Position, PG, remove=FALSE)
# rep 4
annot_4 = read_excel("201201_biofilm_annotation_australian_strains.xlsx", sheet = 'rep4') %>% 
  rename(PG = Plate) %>% 
  mutate(Strain = case_when(Strain == 'OP50' ~ '0',
                            TRUE ~ Strain),
         Strain = as.numeric(Strain)) %>%
  unite('ID', Position, PG, remove=FALSE)


## let's see which category we chose 
annot = annot_1 %>% bind_rows(annot_2,annot_3,annot_4) 

# there are some "lawn_finished" that I'm going to change as to 'normal'
annot %>% count(Annotation_50mM)


annot = annot %>% 
  mutate(Annotation_0mM = case_when(Annotation_0mM == 'lawn_finished' ~ 'normal',
                                    TRUE ~ Annotation_0mM),
         Annotation_50mM = case_when(Annotation_50mM =='lawn_finished' ~ 'normal',
                                     TRUE ~ Annotation_50mM)) %>%
  filter(Annotation_0mM != 'n/a') %>% 
  mutate(pheno_0mM = case_when(Annotation_0mM == 'normal' ~ 0,
                               Annotation_0mM == 'biofilm' ~ 1,
                               Annotation_0mM == 'super_biofilm'~2),
         pheno_50mM = case_when(Annotation_50mM == 'normal' ~ 0,
                               Annotation_50mM == 'biofilm' ~ 1,
                               Annotation_50mM == 'super_biofilm'~2))


# calculate mean and median
annot.sum = annot %>% group_by(ID, PG, Well) %>% 
  summarise(Mean_0 = mean(pheno_0mM, na.rm = T),
            Median_0 = median(pheno_0mM, na.rm = T),
            Mean_50  = mean(pheno_50mM, na.rm = T),
            Median_50 = median(pheno_50mM, na.rm = T),
            Med_0_round= round(Median_0,0),
            Med_50_round = round(Median_50,0))

annot.sum = annot.sum %>% mutate(Annotation_0mM = case_when(Med_0_round == 0 ~ 'normal',
                                                Med_0_round == 1 ~ 'biofilm',
                                                Med_0_round == 2 ~ 'super_biofilm'),
                     Annotation_50mM = case_when(Med_50_round == 0 ~ 'normal',
                                                Med_50_round == 1 ~ 'biofilm',
                                                Med_50_round == 2 ~ 'super_biofilm'))


# THIS METADATA HAS THE MEDIAN OF THE PHENOTYPE ROUNDED UP
biofilm_strain_annotation = biofilm_strain_annotation %>% select(-Annotation_0mM, -Annotation_50mM) %>% 
  left_join(annot.sum %>% select(ID, Annotation_0mM, Annotation_50mM)) %>% 
  select(ID:Strain, Annotation_0mM, Annotation_50mM, notes)


# save the metadata
write.xlsx(biofilm_strain_annotation, '201201_biofilm_annotation_australian_strains_AVERAGE.xlsx')



# join everything to have names, IDs and all that stuff in one table
data_filt_full = data_filt_full %>%
  left_join(biofilm_strain_annotation,by=c('ID','PG')) %>%
  mutate(Metf = factor(Metf, levels = c('0mM', '50mM'))) %>% 
  select(PG, Well, Metf, Annotation_0mM , Annotation_50mM,everything()) %>% 
  replace_na(list('notes' = 'normal'))
# 
# data_filt_full %>% filter(ID=='13_PG2') %>% select(notes, everything()) %>% view()
# 
# data_filt_full %>% filter(notes=='contaminated') %>% select(notes, everything()) %>% view()

# exploratory plots -------------------------------------------------------



# test
# Strain 0 is OP50
data_filt_full %>%
  filter(Strain == 0) %>% 
  ggplot(aes(x = Metf, y = Mean_Intensity, fill = Metf)) +
  geom_boxplot() +
  stat_summary_bin(fun.data = "mean_cl_boot", colour = "red", size = 0.4) +
  geom_point(aes(color = Replicate), position = position_jitterdodge()) +
  facet_wrap(~PG)

ggsave(here('exploration','AUS_control_boxplots.pdf'),width = 10,height = 8)

# test
# Strain 0 is OP50
data_filt_full %>%
  # filter(Strain == 0) %>% 
  ggplot(aes(x = Metf, y = Mean_Intensity, fill = Metf)) +
  geom_boxplot() +
  stat_summary_bin(fun.data = "mean_cl_boot", colour = "red", size = 0.4) +
  # geom_point(aes(color = Replicate), position = position_jitterdodge()) +
  facet_wrap(~PG)


test_0 = data_filt_full %>% 
  unite(test_PG, PG, Metf, remove = F) %>% 
  filter(Metf == '0mM') %>% 
  drop_na(Mean_Intensity,test_PG)

test_0 %>% ggplot(aes(x = test_PG,y=Mean_Intensity, fill = test_PG)) + geom_boxplot()


test_50 = data_filt_full %>% 
  unite(test_PG, PG, Metf, remove = F) %>% 
  filter(Metf == '50mM') %>% 
  drop_na(Mean_Intensity,test_PG) %>% 
  mutate(test_PG = as.factor(test_PG))

test_50 %>% ggplot(aes(x = test_PG,y=Mean_Intensity, fill = test_PG)) + geom_boxplot()



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

data_filt = data_filt_full

data_filt_outs = data_filt %>% 
  group_by(PG,Well,Metf,ID,Annotation_0mM,Annotation_50mM) %>% 
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
write.csv(data_filt_mad, here('exploration', 'AUS_raw_data_filtered.csv'), quote = F, row.names = F)

# save statistics list
list_of_datasets = list('data_filtered' = data_filt_mad)
write.xlsx(list_of_datasets, here('exploration', 'AUS_raw_data_filtered.xlsx'), colNames = T, rowNames = F) 


# Exploration plots -------------------------------------------------------

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




# well order from 1 to 12, with package naturalsort
lvls = naturalsort::naturalsort(unique(data_filt$Well))

# plot number of worms per plate and replicate
data_filt_mad %>%
  # filter(Probability_of_SingleWorm > 0.65) %>%
  filter(PG == 'PG1') %>%
  mutate(Well = factor(Well, levels = lvls)) %>%
  group_by(Well, PG, Replicate, Metf) %>%
  summarise(N = n()) %>%
  ggplot(aes(x = Replicate, y = N, colour = Metf, fill = Metf)) +
  # geom_point(position = position_jitterdodge(), color = 'black') +
  geom_bar(stat = 'identity', width=.5, position = "dodge") +
  geom_hline(yintercept = 3, linetype = 'dashed', colour = 'grey50') +
  facet_wrap(~Well*PG, ncol = 12) +
  theme_classic()


# iterate and save number of worms per plate, well, condition, and replicate
plates = unique(data_filt_mad$PG)

for (plate in plates){
  p = data_filt_mad %>%
    # filter(Probability_of_SingleWorm > 0.65) %>%
    filter(PG == plate) %>%
    mutate(Well = factor(Well, levels = lvls)) %>%
    group_by(ID, Well, PG, Replicate, Metf) %>%
    summarise(N = n()) %>%
    ggplot(aes(x = Replicate, y = N, fill = Metf)) +
    # geom_point(position = position_jitterdodge(), color = 'black') +
    geom_bar(stat = 'identity', width=.85, position = "dodge", color = 'black') +
    geom_hline(yintercept = 3, linetype = 'dashed', colour = 'grey50') +
    scale_fill_manual(values = c('#F5B607', '#0B97D4')) +
    facet_wrap(~ID, ncol = 12) +
    theme_classic()
  ggsave(file = here('exploration/Number_of_worms', paste0('N_worms_',plate,'_replicates.pdf')),
         width = 11, height = 7)
}



# number of worms per condition and plate (reps joined)
for (plate in plates){
  p = data_filt_mad %>%
    # filter(Probability_of_SingleWorm > 0.6) %>%
    filter(PG == plate) %>%
    mutate(Well = factor(Well, levels = lvls)) %>%
    group_by(ID,Well, PG, Metf) %>%
    summarise(N = n()) %>%
    ggplot(aes(x = Metf, y = N, fill = Metf)) +
    # geom_point(position = position_jitterdodge()) +
    geom_bar(stat = 'identity', width=.85, position = "dodge", color = 'black') +
    # ylim(0,25) +
    geom_hline(yintercept = 3, linetype = 'dashed', colour = 'grey50') +
    scale_fill_manual(values = c('#F5B607', '#0B97D4')) +
    facet_wrap(vars(ID), ncol = 12) +
    theme_classic()
  ggsave(file = here('exploration/Number_of_worms', paste0('N_worms_',plate,'_condition.pdf')),
         width = 11, height = 9)
}




### EXPLORATION ####

# # # # # # # # # # # #
### Mean Intensity ####
# # # # # # # # # # # #

IDs = c('17_PG1')
IDs = c('57_PG1','16_PG1','7_PG1','8_PG1')
# plot Mean intensity by itself (it looks good!)
data_filt_mad %>% 
  # filter(Probability_of_SingleWorm > 0.7) %>%
  # filter(Strainname %in% strains) %>%
  filter(ID %in% IDs) %>% 
  ggplot(aes(y = Mean_Intensity , x = ID, fill = Metf)) +
  geom_boxplot(outlier.colour = 'red') +
  geom_point(position = position_jitterdodge(jitter.width = 0.3)) +
  scale_fill_manual(values = c('#F5B607', '#0B97D4')) +
  theme_light() 

# ggsave(file = here('exploration', 'Mean_Intensity_test.pdf'),
#        width = 11, height = 7)



## plot every strain in a single page

strains = data_filt_mad %>%
  select(ID) %>%
  drop_na(ID) %>%
  t %>% as.character %>% unique

plot_list = list()
for (i in 1:length(strains)){
  p = data_filt_mad %>% 
    filter(ID == strains[i]) %>%
    ggplot(aes(y = Mean_Intensity, x = ID, fill = Metf, group = Metf)) +
    geom_boxplot(outlier.colour = 'pink') +
    geom_point(aes(color = Replicate, size = Size_in_pixels), 
               position = position_jitterdodge(jitter.width = 0.3)) +
    # facet_wrap(~PG) +
    scale_fill_manual(values = c('#F5B607', '#0B97D4')) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 0, hjust  = 1))
  
  plot_list[[i]] = p
}

# save multiple plots into one big file
ggsave(file = here('analysis', 'Mean_Intensity_ALL_STRAINS.pdf'), 
       marrangeGrob(grobs = plot_list, nrow = 1, ncol = 1),
       width = 11, height = 7)




# # # # # # # # # # # # # #
### Summary statistics ####
# # # # # # # # # # # # # #


# this piece of code calculates FC per replicate, and then calculate the mean 
# of the FCs, its SD and also the median
FC_means = data_filt_mad %>% 
  group_by(PG, Replicate, Well, Metf, ID, Annotation_0mM, Annotation_50mM) %>% 
  summarise(Mean = mean(Mean_Intensity)) %>% 
  pivot_wider(names_from = Metf, values_from = Mean, names_prefix = 'Metf_') %>% 
  mutate(FC = Metf_50mM/Metf_0mM) %>%
  group_by(PG, Well, ID, Annotation_0mM, Annotation_50mM)  %>% 
  summarise(Mean_FC = mean(FC, na.rm = TRUE),
            SD_FC = sd(FC, na.rm = TRUE),
            Median_FC = median(FC, na.rm = TRUE))

FC_means %>% 
  ungroup %>% 
  arrange(desc(Mean_FC)) %>% 
  mutate(ID = factor(ID, levels = ID)) %>% 
  ggplot(aes(x = ID, y = Mean_FC)) +
  geom_pointrange(aes(ymin = Mean_FC - SD_FC, ymax = Mean_FC + SD_FC)) +
  geom_hline(yintercept = 1) +
  theme_classic() +
  labs(x='Strains',
       y='FC of 50mM/0mM (Mean +- SD)') +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size=4))

ggsave(here('exploration','FC_ordered.pdf'),height = 7,width = 17)


# save a list of unique IDs 
FC_means_unique = FC_means %>% 
  ungroup %>% 
  distinct(ID, .keep_all = T)

write_csv(FC_means, here('analysis', 'FC_means.csv'))

write_csv(FC_means_unique, here('analysis', 'FC_means_unique.csv'))



sum.stats = data_filt_mad %>% 
  group_by(PG, Metf, ID, Annotation_0mM, Annotation_50mM) %>% 
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

library(broom)

# some samples are bad, so remove them for now
removals = c('7_PG1','8_PG1','7_PG2','28_PG3','52_PG3')

# per plate and ID
stat_res = data_filt_mad %>%
  filter(!(ID %in% removals)) %>% 
  group_by(PG, ID, Annotation_0mM, Annotation_50mM) %>%
  nest() %>%
  mutate(models = map(data, lm, formula = 'Mean_Intensity ~ Metf')) %>%
  mutate(results = map(.f = tidy, .x = models)) %>%
  select(PG, Annotation_0mM, Annotation_50mM, results) %>%
  unnest(cols = c(results)) %>%
  filter(term != '(Intercept)') %>% 
  ungroup %>% # ungroup before passing the FDR correction or it will correct line by line
  mutate(FDR = p.adjust(p.value, method = 'fdr'),
         FDR_stars = gtools::stars.pval(FDR))

stat_res = data_filt_mad %>%
  group_by(ID, PG, Metf, Well) %>%
  summarise(Mean = mean(Mean_Intensity, na.rm = TRUE)) %>%
  pivot_wider(names_from = c(Metf), values_from = Mean, names_prefix = 'Met_') %>%
  mutate(FC = Met_50mM/Met_0mM,
         log2FC = log2(FC)) %>%
  select(ID, PG, Well, FC, log2FC, Met_0mM, Met_50mM) %>%
  left_join(stat_res) %>%
  ungroup %>%
  select(PG, Well, ID, Annotation_0mM, Annotation_50mM, FC, log2FC,Met_0mM, Met_50mM, estimate:FDR_stars)


# metadata for stats
tables = c('Stats_per_strain')
description = c('Statistics with individual duplicates')
stats_metadata = data.frame(tables, description)

# save statistics list
list_of_datasets = list('metadata' = stats_metadata, 'Stats_per_plate' = stat_res)

write.xlsx(list_of_datasets, here('analysis', 'worm_imaging_stats.xlsx'), colNames = T, rowNames = F) 
write.xlsx(list_of_datasets, here('exploration', 'worm_imaging_stats.xlsx'), colNames = T, rowNames = F) 



# mean plots
stat_res %>% 
  mutate(ID = fct_reorder(ID, desc(Met_0mM))) %>% 
  drop_na(Annotation_0mM) %>% 
  ggplot(aes(x=ID, y = Met_0mM, colour = Annotation_0mM)) +
  geom_point()+
  theme_classic() +
  labs(x='strains',
       y = 'Mean Pacs2::GFP activation at 0mM metformin') +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size=4))

ggsave(here('exploration','Met0mM_ordered.pdf'),height = 7,width = 17)


stat_res %>% 
  drop_na(Annotation_50mM) %>% 
  mutate(ID = fct_reorder(ID, desc(Met_50mM))) %>% 
  ggplot(aes(x=ID, y = Met_50mM, colour = Annotation_50mM)) +
  geom_point()+
  theme_classic() +
  labs(x='strains',
       y = 'Mean Pacs2::GFP activation at 50mM metformin') +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size=4))

ggsave(here('exploration','Met50mM_ordered.pdf'),height = 7,width = 17)













