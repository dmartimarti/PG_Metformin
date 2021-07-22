
# libraries ---------------------------------------------------------------
library(tidyverse)
library(readxl)
library(broom)
library(openxlsx)
library(here)

theme_set(theme_classic())

# Read data ---------------------------------------------------------------


data = read_csv('Output/Summary.csv', quote = "\"") %>%
  rename(AUC_raw = `595nm_f_AUC`) %>% # `750nm_f_logAUC` data column is what we need for logAUC values
  mutate(Strain = as.character(Strain),
         Row = str_match_all(Well,'[:digit:]{1,}'), #Get plate row names from well name. 
         Col = str_match_all(Well,'[:alpha:]{1,}'), #Get plate column names from well name
         Row = factor(Row, levels = 1:12), #Make them categorical variable with set order them
         Col = factor(Col, levels = LETTERS[1:8]),
         Strain = as.factor(Strain)) %>% #Change Type column coding
  drop_na(Strain) %>% 
  unite(ID, Strain, Plate, Well, remove = F) %>%
  select(-Replicate) %>% 
  select(Strain, ID, Replicate = Bio_rep, Metformin_mM = Met, Plate, Well, Row, Col, AUC_raw)


metadata = read_excel("D:/MRC_Postdoc/Pangenomic/pangenome_analysis/metadata/MAIN_metadata.xlsx", 
                            sheet = "metadata")



### IMPORTANT: After running this chunck, I didn't detect any major outliers, so 
#  I decided to not use this outlier detection tool 

# outlier detection -------------------------------------------------------
# 
# the process of outlier detection has these steps:
#   1. use DBSCAN with a wide table (select optimal eps with kNNdistplot)
#   2. check and get outliers from analysis for further processing and outlier refining
#   3. calculate z-score and coefficient of variation (CV)
#   4. use the samples that deviate greatly from a threshold (0.8)
#   5. to select the replicate that is the outlier, calculate pairwise differences between their z-scores
#   6. get the differences that are larger than 1.3 (see plot to check and select threshold)
#   7. select only those samples that had only 1 bad replicate 
#   8. remove bad replicates and clean the dataset
# 
# clean up of data
data.wide = data %>%
  select(ID, Metformin_mM, Replicate, AUC_raw) %>%
  pivot_wider(names_from = Replicate, values_from = AUC_raw, names_prefix = 'rep_')


library(plotly)


plot_ly(x=data.wide$rep_1, y=data.wide$rep_2, z=data.wide$rep_3,
        type="scatter3d", mode="markers",
        color = as.factor(data.wide$Metformin_mM))



# #   1. use DBSCAN with a wide table (select optimal eps with kNNdistplot)
# 
# library(fpc)
# library(dbscan)
# 
# # calculate the optimal eps for db
# dbscan::kNNdistplot(data.wide[,3:5], k =  5)
# abline(h = 1.4, lty = 2)
# 
# db = fpc::dbscan(data.wide[,3:5], eps = 2, MinPts = 5,
#                  method = 'hybrid')
# 
# db
# 
# #   2. check and get outliers from analysis for further processing and outlier refining
# # plot outliers
# plot_ly(x=data.wide$rep_1, y=data.wide$rep_2, z=data.wide$rep_3, 
#         type="scatter3d", mode="markers",
#         color = as.factor(db$cluster))
# 
# 
# outliers = db$cluster
# 
# out_IDs = data.wide %>% cbind(outliers) %>% 
#   as_tibble %>% 
#   filter(outliers == 0) %>% 
#   separate(ID, into = c('ID', 'Plate', 'Well'), sep = '_') %>% 
#   pull(ID)
# 
# 
# out_dup = data %>% 
#   filter(Strain %in% out_IDs) %>% 
#   group_by(Strain, Metformin_mM) %>% 
#   mutate(zscore = (AUC_raw-mean(AUC_raw))/sd(AUC_raw)) %>% 
#   ungroup 
# 
# 
# #   3. calculate z-score and coefficient of variation (CV)
# # refine selection by calculating the coefficient of variation
# cv = out_dup %>% 
#   group_by(Strain, Metformin_mM, ID) %>% 
#   # mutate(AUC_raw = log2(AUC_raw)) %>% 
#   summarise(Mean = mean(AUC_raw),
#             SD = sd(AUC_raw)) %>% 
#   mutate(CV = SD / Mean)
# 
# #   4. use the samples that deviate greatly from a threshold (0.8)
# # mark a global threshold of 0.8
# cv %>% 
#   ggplot(aes(x = fct_reorder(ID, CV), y = CV)) +
#   geom_point() +
#   geom_hline(yintercept=0.8) +
#   facet_wrap(~Metformin_mM) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
# 
# ggsave(here('summary', 'outliers_coef_variation.pdf'), height = 9, width = 10)
# 
# # check data to be sure the threshold is ok
# out_dup %>% left_join(cv) %>% view
# 
# # get true outliers
# true_outs = cv %>% ungroup %>%  filter(CV > 0.8) %>% 
#   select(Strain, Metformin_mM, ID) %>% 
#   unite(new_ID, ID, Metformin_mM, remove=F) %>% pull(new_ID)
# 
# # filter the previous table
# out_dup = out_dup %>% 
#   unite(new_ID, ID, Metformin_mM, remove=F) %>% 
#   filter(new_ID %in% true_outs) %>% 
#   select(-new_ID) %>% 
#   left_join(cv)
# 
# #   5. to select the replicate that is the outlier, calculate pairwise differences between their z-scores
# out_calc = out_dup %>% 
#   select(ID, Strain, Metformin_mM, Replicate, zscore) %>% 
#   pivot_wider(names_from = c('Replicate'), names_prefix = 'rep_', values_from = zscore) %>% 
#   mutate(dif_1_2 = abs(rep_1 - rep_2),
#          dif_2_3 = abs(rep_2 - rep_3),
#          dif_1_3 = abs(rep_1 - rep_3)) %>% 
#   select(ID, Strain, Metformin_mM, dif_1_2:dif_1_3) %>% 
#   pivot_longer(dif_1_2:dif_1_3, names_to = 'Groups', values_to = 'difference') %>% 
#   mutate(Replicate = case_when(Groups == 'dif_1_2' ~ 3,
#                                Groups == 'dif_2_3' ~ 1,
#                                Groups == 'dif_1_3' ~ 2)) 
# #   6. get the differences that are larger than 1.3 (see plot to check and select threshold)
# out_calc %>% 
#   ggplot(aes(x = fct_reorder(ID, difference), y = difference)) +
#   geom_point() + facet_wrap(~Metformin_mM) +
#   geom_hline(yintercept = 1.2)
# 
# #   7. select only those samples that had only 1 bad replicate 
# IDs_bad = out_calc %>% filter(difference < 1.2) %>% 
#   count(ID,Metformin_mM) %>% arrange(desc(n)) %>% 
#   filter(n == 1) %>%  ## keep only the ones that are 
#   pull(ID)
# 
# #   8. remove bad replicates and clean the dataset
# ### THIS IS THE IMPORTANT LIST
# # this list has the IDs of the bad replicates that we need to remove from the analysis
# outliers_final = out_calc %>% filter(difference < 1.2) %>% 
#   filter(ID %in% IDs_bad) %>% 
#   select(ID, Strain, Metformin_mM, Replicate)
# 
# 
# # CLEAN THE DATA FROM THE OUTLIERS FOUND
# data.clean = data %>% anti_join(outliers_final)
# 
# ## BONUS: check that the distribution is better
# # clean up of data
# data.wide.clean = data.clean %>% 
#   select(ID, Metformin_mM, Replicate, AUC_raw) %>% 
#   pivot_wider(names_from = Replicate, values_from = AUC_raw, names_prefix = 'rep_')
# 
# # check that it looks better than before
# plot_ly(x=data.wide.clean$rep_1, y=data.wide.clean$rep_2, z=data.wide.clean$rep_3, 
#         type="scatter3d", mode="markers",
#         color = as.factor(data.wide.clean$Metformin_mM))
# 


data.clean = data



# Boxplots -------------------------------------------------------

plates = c(1,2,3,4)

for (plate in plates){
  data.clean %>% 
    mutate(Metformin_mM = factor(Metformin_mM, levels = c(0,50,100,200))) %>% 
    filter(Plate == plate) %>% 
    ggplot(aes(x = Metformin_mM, y = AUC_raw, fill = Metformin_mM)) +
    geom_boxplot() +
    facet_wrap(~Strain, ncol = 12)
  print(paste0('printing plate ', as.character(plate)))
  # save multiple plots into one big file
  ggsave(file = here('summary', paste0('plate_',as.character(plate),'_boxplot.pdf')),
         width = 13, height = 8)
}

data.clean %>% filter(Plate == 2) %>% view()



# cummulative AUCs --------------------------------------------------------



# long transformation of data
data.sum = data.clean %>% 
  # anti_join(bad) %>% 
  # generate normalised AUCs dividing each value by the 0 metf
  pivot_wider(names_from = c(Replicate, Metformin_mM), values_from = AUC_raw, names_prefix = 'rep_') %>% 
  mutate(AUC_1_50 =  rep_1_50/ rep_1_0,
         AUC_1_100 = rep_1_100/rep_1_0,
         AUC_1_200 = rep_1_200/rep_1_0,
         AUC_3_50 =  rep_3_50/ rep_3_0,
         AUC_2_50 =  rep_2_50/ rep_2_0,
         AUC_2_100 = rep_2_100/rep_2_0,
         AUC_2_200 = rep_2_200/rep_2_0,
         AUC_3_100 = rep_3_100/rep_3_0,
         AUC_3_200 = rep_3_200/rep_3_0
  ) %>% 
  # uncsum = unweighted cummulative sum
  # wcsum = weighted cummulative sum
  select(Strain:Col, AUC_1_50:AUC_3_200) %>% 
  mutate(uncsum_1 = AUC_1_50 + AUC_1_100 + AUC_1_200,
         uncsum_2 = AUC_2_50 + AUC_2_100 + AUC_2_200,
         uncsum_3 = AUC_3_50 + AUC_3_100 + AUC_3_200,
         wcsum_1 = AUC_1_50 + (AUC_1_100)*2 + (AUC_1_200)*4,
         wcsum_2 = AUC_2_50 + (AUC_2_100)*2 + (AUC_2_200)*4,
         wcsum_3 = AUC_3_50 + (AUC_3_100)*2 + (AUC_3_200)*4) %>% 
  select(ID, uncsum_1:wcsum_3) %>%
  pivot_longer(-ID, names_to = 'Measure', values_to = 'AUC') %>% 
  separate(col = Measure, sep = '_', into = c('Measure', 'Replicate')) %>% 
  # filter(Replicate != 1) %>%
  group_by(ID, Measure) %>% 
  summarise(Mean = mean(AUC, na.rm = T),
            SD = sd(AUC, na.rm = T)) %>% 
  ungroup




data.sum %>% 
  filter(Mean < 10) %>% 
  group_by(Measure) %>% 
  arrange(Mean) %>%
  mutate(ID = factor(ID, levels = ID)) %>% 
  ggplot(aes(x = ID, y = Mean)) +
  # geom_point() +
  geom_pointrange(aes(ymin = Mean - SD, ymax = Mean + SD)) +
  facet_wrap(~Measure, ncol = 1)


sum.stats = data.clean %>% 
  filter(Plate %in% c(1,2,3,4,5,6)) %>% 
  group_by(Strain, ID, Plate, Well, Metformin_mM) %>% 
  summarise(Mean = mean(AUC_raw, na.rm = T),
            SD = sd(AUC_raw, na.rm = T)) %>% 
  pivot_wider(names_from = Metformin_mM, names_prefix = 'Bact_metf_', values_from = c(Mean, SD))



# filter out the exp strains
data.sum = data.sum %>%
  separate(ID, into = c('Strain', 'Plate', 'Well'), sep = '_') 

data.sum = data.sum %>% 
  mutate(Plate = as.integer(Plate)) %>% 
  left_join(sum.stats) %>% 
  rename(Bact_score_mean = Mean,
         Bact_score_sd = SD) %>% 
  select(Strain, ID, Plate, Well, Measure, everything(), -ID)


data.sum = data.sum %>% left_join(metadata %>% filter(Origin == 'AUS') %>% 
  select(Strain = ID, phylogroup))

# fill empty fields with non ecoli labels
data.sum = data.sum %>% 
  mutate(phylogroup = case_when(is.na(phylogroup) ~ 'Non Escherichia',
                                TRUE ~ phylogroup)) 

# save info
list_of_datasets = list('unweighted' = data.sum %>% filter(Measure == 'uncsum'), 
                        'weighted' = data.sum %>% filter(Measure == 'wcsum'))

write.xlsx(list_of_datasets, 'Growth_resistance_summaryStats_AUS.xlsx', colNames = T, rowNames = T) 


data.sum %>% 
  # filter(Measure == 'uncsum') %>% 
  filter(phylogroup != 'Non Escherichia') %>% 
  mutate(Measure = case_when(Measure == 'uncsum' ~ 'Unweighted',
                             TRUE ~ 'Weighted')) %>%
  ggplot(aes(y = Bact_score_mean, x = fct_reorder(Strain, Bact_score_mean), color = Bact_score_mean)) +
  geom_point() +
  facet_wrap(~Measure, nrow=2, scales = 'free_y')
  
ggsave(here('summary', 'growth_score.pdf'), height = 8, width = 10)

data.sum %>% 
  # filter(Measure == 'uncsum') %>% 
  filter(phylogroup != 'Non Escherichia') %>% 
  mutate(Measure = case_when(Measure == 'uncsum' ~ 'Unweighted',
                             TRUE ~ 'Weighted')) %>% 
  select(Strain, Measure, Bact_score_mean) %>% 
  pivot_wider(names_from = 'Measure', values_from = Bact_score_mean) %>% 
  ggplot(aes(x = Unweighted, y = Weighted)) +
  geom_smooth() +
  geom_point()

# Growth curves (inherited from Jen) --------------------------------------


### TIME SERIES, growth assays were performed in 60-well format

# Get timeseries data
time_data = read_csv('Output/Timeseries.csv', quote = "\"") %>%
  filter(Data == '595nm_f') %>%
  select(-Pattern, -Reader) %>%
  gather(Time_s, OD, matches('\\d')) %>%
  filter(!is.na(OD)) %>% # Remove empty values if there are missmatches
  filter(!is.na(Strain)) %>% # Remove empty wells without strains
  mutate(Time_s = as.numeric(Time_s),
         Time_h = Time_s/3600,
         Row = str_match_all(Well,'[:digit:]{1,}'), #Get plate row names from well name.
         Col = str_match_all(Well,'[:alpha:]{1,}'), #Get plate column names from well name
         Row = factor(Row, levels = 1:12), #Make them categorical variable with set order them
         Col = factor(Col, levels = LETTERS[1:8]),
         Strain = as.factor(Strain)) %>%
  select(-File, -Data, -Media) %>%
  rename(Metformin_mM	= Met)

# get factor sort order depending on the wells
lvls = naturalsort::naturalsort(unique(time_data$Well))

# Divde per plate, since it's easier to plot them separate

AUS1 = time_data %>%
  filter(Plate == '1')

AUS2 = time_data %>%
  filter(Plate == '2')

AUS3 = time_data %>%
  filter(Plate == '3')

AUS4 = time_data %>%
  filter(Plate == '4')


## Check if your biological replicates have the same timing (Hours), sometimes it is one second off despite using the same Biospa session
times = AUS2 %>% filter(Strain == '66' & Bio_rep == '3' & Metformin_mM == '0') %>% select(Time_h) %>% t %>% as.character
timescheckrep2 = AUS2 %>% filter(Strain == '66' & Bio_rep == '2' & Metformin_mM == '0') %>% select(Time_h) %>% t %>% as.character
timescheckrep3 = AUS2 %>% filter(Strain == '66' & Bio_rep == '1' & Metformin_mM == '0') %>% select(Time_h) %>% t %>% as.character

times = AUS1 %>% filter(Strain == '2' & Bio_rep == '2' & Metformin_mM == '0') %>% select(Time_h) %>% t %>% as.character
timescheckrep2 = AUS1 %>% filter(Strain == '2' & Bio_rep == '1' & Metformin_mM == '0') %>% select(Time_h) %>% t %>% as.character
timescheckrep3 = AUS1 %>% filter(Strain == '2' & Bio_rep == '3' & Metformin_mM == '0') %>% select(Time_h) %>% t %>% as.character

#times = AUS3 %>% filter(Strain == '127' & Bio_rep == '1' & Metformin_mM == '0') %>% select(Time_h) %>% t %>% as.character
#timescheckrep2 = AUS3 %>% filter(Strain == '127' & Bio_rep == '2' & Metformin_mM == '0') %>% select(Time_h) %>% t %>% as.character
#timescheckrep3 = AUS3 %>% filter(Strain == '127' & Bio_rep == '3' & Metformin_mM == '0') %>% select(Time_h) %>% t %>% as.character

#times = AUS4 %>% filter(Strain == '187' & Bio_rep == '1' & Metformin_mM == '0') %>% select(Time_h) %>% t %>% as.character
#timescheckrep2 = AUS4 %>% filter(Strain == '187' & Bio_rep == '2' & Metformin_mM == '0') %>% select(Time_h) %>% t %>% as.character
#timescheckrep3 = AUS4 %>% filter(Strain == '187' & Bio_rep == '3' & Metformin_mM == '0') %>% select(Time_h) %>% t %>% as.character


# correct the timing if needed
AUS1 = AUS1 %>%
  group_by(Strain, Bio_rep, Metformin_mM) %>%
  mutate(Time_h = times) %>%
  ungroup

AUS2 = AUS2 %>%
  group_by(Strain, Bio_rep, Metformin_mM) %>%
  mutate(Time_h = times) %>%
  ungroup

# summarise data for bio replicates, use either original data tibble or the time-fixed tibble

tsum = AUS4 %>%
  group_by(Strain, Metformin_mM, Well, Time_h) %>%
  summarise(Mean = mean(OD), SD = sd(OD), SE = SD/sqrt(length(OD))) %>%
  ungroup

tsum$Strain <- sub("^", "str", tsum$Strain )

####### To plot the curves per plate

tsum %>%
  mutate(Well = factor(Well, levels = lvls), Metformin_mM = as.factor(Metformin_mM), Time_h = as.numeric(Time_h)) %>%
  ggplot( aes(x = Time_h, y = Mean, fill = Metformin_mM, color = Metformin_mM)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.2) +
  geom_line(size = 1) +
  scale_x_continuous(breaks = seq(0, 24, by = 6)) +
  labs(x = 'Time, h',
       y = 'O.D.') +
  labs(fill = "Metformin_mM") +
  facet_wrap(vars(Strain), ncol = 10) +
  scale_colour_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey50"))

quartz.save(file = 'AUS_P4.pdf',
            type = 'pdf', dpi = 300, height = 6, width = 10)









