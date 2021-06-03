
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


# outlier detection -------------------------------------------------------

library(OutlierDetection)

# filter out the IDs (unique strains) that have a growth less than 0.05 at metf == 0
bad = data %>% 
  filter(AUC_raw < 0.05, Metformin_mM == 0) %>% 
  arrange(ID) %>% distinct(ID) %>% t %>% as.character()

# clean up of data
data.wide = data %>% 
  # DATA IMPUTATION by minimum non-zero value
  # mutate(AUC_raw = case_when(AUC_raw == 0 ~ 0.0001,
  #                            AUC_raw != 0 ~ AUC_raw)) %>% 
  # filter(!(ID %in% bad)) %>% 
  select(ID, Metformin_mM, Replicate, AUC_raw) %>% 
  pivot_wider(names_from = Replicate, values_from = AUC_raw, names_prefix = 'rep_')

# met0_filt = data %>% filter(Metformin_mM == 0) %>%
#   select(ID, Replicate, AUC_raw) %>%
#   pivot_wider(names_from = Replicate, values_from = AUC_raw, names_prefix = 'rep_') 

# get the list of weird strains
# weird_st = met0_filt %>% filter(rep_1 <= thr | rep_2 <= thr | rep_3 <= thr) %>%
#   select(ID) %>% as.vector %>% t %>% as.character

# remove names
# met0_filt = met0_filt %>% filter(!ID %in% weird_st)
X = data.wide[,3:5]

# master method
outs = OutlierDetection(X, k = 5, cutoff = 0.98, Method = "euclidean", 
                        rnames = data.wide[,1], dispersion = TRUE)

outs

# how many outlier strains have duplicates
outliers = data.wide[outs$`Location of Outlier`,] %>% 
  mutate(ID = str_sub(ID, 1,7)) %>% 
  distinct(ID) %>% 
  t %>% as.character()

# 64 strains with duplicates
sum(outliers %in% dups)

# filter duplicated strains in the outlier group
out_dup = data %>% 
  filter(Strain %in% outliers[outliers %in% dups]) %>% 
  group_by(Strain, Metformin_mM) %>% 
  mutate(zscore = (AUC_raw-mean(AUC_raw))/sd(AUC_raw)) %>% 
  ungroup 

# use the outliers function from spatialEco to calculate 
# the modified z-score, that is better suited to find 
# outliers 

# initialise thresholds
thr1 = 9
thr2 = 9
thr3 = 7
thr4 = 3

out_dup = out_dup %>% 
  group_by(Strain, Metformin_mM) %>% 
  mutate(out = spatialEco::outliers(AUC_raw)) %>% 
  ungroup %>% 
  mutate(pass = case_when(Metformin_mM == 0 & abs(out) > thr1 ~ 'OUT',
                          Metformin_mM == 0 & abs(out) <= thr1 ~ 'pass',
                          Metformin_mM == 50 & abs(out) > thr2 ~ 'OUT',
                          Metformin_mM == 50 & abs(out) <= thr2 ~ 'pass',
                          Metformin_mM == 100 & abs(out) > thr3 ~ 'OUT',
                          Metformin_mM == 100 & abs(out) <= thr3 ~ 'pass',
                          Metformin_mM == 200 & abs(out) > thr4 ~ 'OUT',
                          Metformin_mM == 200 & abs(out) <= thr4 ~ 'pass')) %>% 
  arrange(Strain, Metformin_mM)

# how many replicates we have per condition?
# at least 4, so far so good, we can remove the outliers detected previously

out_dup %>% 
  filter(pass != 'OUT') %>% 
  group_by(Strain, Metformin_mM) %>% 
  count() %>% arrange(n)


### do the same, but with non-duplicates
out_single = data %>% 
  filter(Strain %in% outliers[!(outliers %in% dups)]) %>% 
  group_by(Strain, Metformin_mM) %>% 
  mutate(zscore = (AUC_raw-mean(AUC_raw))/sd(AUC_raw)) %>% 
  ungroup


out_single = out_single %>% 
  group_by(Strain, Metformin_mM) %>% 
  mutate(out = spatialEco::outliers(AUC_raw)) %>% 
  ungroup %>% 
  mutate(pass = case_when(Metformin_mM == 0 & abs(out) > thr1 ~ 'OUT',
                          Metformin_mM == 0 & abs(out) <= thr1 ~ 'pass',
                          Metformin_mM == 50 & abs(out) > thr2 ~ 'OUT',
                          Metformin_mM == 50 & abs(out) <= thr2 ~ 'pass',
                          Metformin_mM == 100 & abs(out) > thr3 ~ 'OUT',
                          Metformin_mM == 100 & abs(out) <= thr3 ~ 'pass',
                          Metformin_mM == 200 & abs(out) > thr4 ~ 'OUT',
                          Metformin_mM == 200 & abs(out) <= thr4 ~ 'pass')) %>% 
  arrange(Strain, Metformin_mM)

# at least two samples in each condition, not bad
out_single %>% 
  filter(pass != 'OUT') %>% 
  group_by(Strain, Metformin_mM) %>% 
  count() %>% arrange(n)



out_final = out_dup %>% filter(pass == 'OUT') %>% 
  bind_rows(out_single %>% filter(pass == 'OUT')) %>% 
  select(-zscore:-pass)


# CLEAN THE DATA FROM THE OUTLIERS FOUND
data.clean = data %>% anti_join(out_final)


# remove PG7 and PG8
data.clean = data.clean %>% filter(Plate %in% c(1,2,3,4,5,6))


# # # 
# OUTLIERS NOT DETECTED YET IN THE DATASET # # #



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


# save info
list_of_datasets = list('unweighted' = data.sum %>% filter(Measure == 'uncsum'), 
                        'weighted' = data.sum %>% filter(Measure == 'wcsum'))

write.xlsx(list_of_datasets, 'Growth_resistance_summaryStats_AUS.xlsx', colNames = T, rowNames = T) 



# comparison --------------------------------------------------------------

# this chunk of code is meant to compare the first version of
# score calculation with the second, where I just calculated first 
# the average of each measure at a certain concentration and then
# I calculated the score

merge = data.sum %>% 
  filter(Measure == 'uncsum') %>% 
  distinct(Strain, .keep_all = T) %>% 
  left_join(nodups) %>% 
  drop_na(Mean, Bact_score_mean) 

merge %>% 
  filter(Bact_score_mean < 5) %>%
  ggplot(aes(x = Mean, y = Bact_score_mean)) +
  geom_smooth(method = 'lm')+
  geom_point() +
  theme_light() +
  labs(x = 'Method 2',
       y = 'Method 1')



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









