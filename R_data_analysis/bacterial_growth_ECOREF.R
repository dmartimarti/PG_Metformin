
# libraries ---------------------------------------------------------------
library(tidyverse)
library(readxl)
library(broom)
library(openxlsx)
library(here)

theme_set(theme_classic())

# Read data ---------------------------------------------------------------


data = read_csv('Summary.csv', quote = "\"") %>%
  rename(AUC_raw = `595nm_f_AUC`) %>% # `750nm_f_logAUC` data column is what we need for logAUC values
  mutate(Strain = as.character(Strain),
         Row = str_match_all(Well,'[:digit:]{1,}'), #Get plate row names from well name. 
         Col = str_match_all(Well,'[:alpha:]{1,}'), #Get plate column names from well name
         Row = factor(Row, levels = 1:12), #Make them categorical variable with set order them
         Col = factor(Col, levels = LETTERS[1:8]),
         Strain = as.factor(Strain)) %>% #Change Type column coding
  unite(ID, Strain, Plate, Well, remove = F) %>%
  select(Strain, ID, Replicate, Metformin_mM, Replicate, Plate, Well, Row, Col, AUC_raw)


# filter empty wells (some of them have growth)
data = data %>% filter(Strain != 'EMPTY')



#### FIX PG6 ####

# according to Jen, 50 and 0 mM in PG6 are switched
# this code fixes it
sub50 = data %>% 
  filter(Plate == 6 & Metformin_mM == 0 & Replicate == 1) %>% 
  mutate(Metformin_mM = 50)
sub50 = data %>% 
  filter(Plate == 6 & Metformin_mM == 50 & Replicate == 1) %>% 
  mutate(Metformin_mM = 0) %>% 
  bind_rows(sub50)

data = data %>% 
  filter(!(Plate == 6 & Replicate == 1 & Metformin_mM %in% c(0,50)))  %>% 
  bind_rows(sub50)

# read strain db database
straindb = read_xlsx('strain_db.xlsx', sheet = 'strain_db')
exp_str = straindb %>% filter(Broadphenotype == 'Evolutionexperiment') %>% select(ID) %>% t %>% as.vector

# filter out the evo strains
data = data %>% filter(!(Strain %in% exp_str))




### duplicates ####

dups = data %>%
  group_by(Strain) %>%
  summarise(N = n()/12) %>%
  arrange(desc(N)) %>%
  drop_na(Strain) %>%
  filter(N > 1) %>% select(Strain) %>% t %>% as.character

num = 5
data %>% 
  filter(Strain == dups[num]) %>% 
  mutate(Metformin_mM = as.factor(Metformin_mM)) %>% 
  ggplot(aes(x = Metformin_mM, y = AUC_raw, group = Metformin_mM, fill = Metformin_mM)) + 
  geom_boxplot() +
  facet_wrap(~Plate)

data %>% 
  filter(Strain == dups[num]) %>% 
  mutate(Metformin_mM = as.factor(Metformin_mM)) %>% 
  ggplot(aes(x = Metformin_mM, y = AUC_raw, group = Metformin_mM, fill = Metformin_mM)) + 
  geom_boxplot() 



# outlier detection -------------------------------------------------------

# the library is no longer active, so this is how you can install it
require(devtools)
# install_version("OutlierDetection", version = "0.1.1", 
#                 repos = "http://cran.us.r-project.org")
# 
# install.packages('https://cran.r-project.org/src/contrib/Archive/OutlierDetection/OutlierDetection_0.1.1.tar.gz',
#                  repos=NULL)

library(OutlierDetection)

# filter out the IDs (unique strains) that have a growth less than 0.05 at metf == 0
bad = data %>% 
  filter(AUC_raw < 0.05, Metformin_mM == 0) %>% 
  filter(!(Strain == 'NT12335')) %>% 
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



# Boxplots -------------------------------------------------------

plates = c(1,2,3,4,5,6)

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



# overlap test ------------------------------------------------------------

## This chunk of code is meant to test how many strains are the same whether 
## we include all reps or only reps2 and 3. 

# first_all = data.sum %>% arrange(Mean) %>% filter(Mean < 15, Measure == 'wcsum') %>%  head(30)
# last_all = data.sum %>%  arrange(Mean) %>% filter(Mean < 15, Measure == 'wcsum') %>%  tail(30)
# 
# first_reps = data.sum %>% arrange(Mean) %>% filter(Mean < 15, Measure == 'wcsum') %>%  head(30)
# last_reps = data.sum %>%  arrange(Mean) %>% filter(Mean < 15, Measure == 'wcsum') %>%  tail(30)
# 
# sum(first_all$ID %in% first_reps$ID)/30
# sum(last_all$ID %in% last_reps$ID)/30

## uncsum
# 67% of top are the same if we include all reps or only reps2 and 3
# 93% of tail are the same 

## wcsum
# 60% of top are the same if we include all reps or only reps2 and 3
# 93% of tail are the same 

## final results: 2 thirds of the top are shared, and by the look of the plot, I don't think the 
## top 40 or top50 would be much different. So, I would include all reps and use the uncsum measure
## to arrange strains by their resistance

# filter out the exp strains
data.sum = data.sum %>%
  separate(ID, into = c('Strain', 'Plate', 'Well'), sep = '_') %>% 
  filter(!(Strain %in% exp_str)) 

data.sum = data.sum %>% 
  mutate(Plate = as.integer(Plate)) %>% 
  left_join(sum.stats) %>% 
  rename(Bact_score_mean = Mean,
         Bact_score_sd = SD) %>% 
  select(Strain, ID, Plate, Well, Measure, everything(), -ID)

# save info
list_of_datasets = list('unweighted' = data.sum %>% filter(Measure == 'uncsum'), 
                        'weighted' = data.sum %>% filter(Measure == 'wcsum'))

write.xlsx(list_of_datasets, 'Growth_resistance_summaryStats.xlsx', colNames = T, rowNames = T) 








# Summary no dups ---------------------------------------------------------


nodups = data.clean %>% 
  # anti_join(bad) %>% 
  filter(!Strain %in% exp_str) %>% 
  filter(!(Strain == 'NT12128' & Replicate == 2 & Metformin_mM == 50)) %>%  # this rep is bad
  # filter(!Strain == 'NT12335') %>%
  # filter(!(Strain == 'NT12335' & Replicate == 1 & Metformin_mM == 100)) %>%
  # filter(!(Strain == 'NT12335' & Replicate == 1 & Metformin_mM == 50)) %>%
  group_by(Strain, Metformin_mM) %>% 
  summarise(Mean = mean(AUC_raw, na.rm = TRUE),
            SD = sd(AUC_raw, na.rm = TRUE),
            Median = median(AUC_raw, na.rm = TRUE),
            SEM = SD/sqrt(n()))

# check how it looks 
nodups %>% 
  filter(Metformin_mM == 200) %>% 
  ungroup %>% 
  arrange(desc(Mean)) %>% 
  mutate(Strain = factor(Strain, levels = Strain)) %>% 
  ggplot(aes(x = Strain, y = Mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM))

# data transformation
nodups = nodups %>%
  select(-c(SD, Median, SEM)) %>% 
  pivot_wider(names_from = Metformin_mM, values_from = Mean, names_prefix = 'Metf_') %>%
  mutate(AUC50 = Metf_50/Metf_0,
         AUC100 = Metf_100/Metf_0,
         AUC200 = Metf_200/Metf_0,
         Mean = AUC50 + AUC100 + AUC200) %>% 
  select(Strain,AUC50,AUC100,AUC200, Mean)

nodups = nodups %>% left_join(sum.stats %>% 
                                ungroup %>% 
                                distinct(Strain, .keep_all = T)) %>% 
  ungroup

# save info
list_of_datasets = list('unweighted' = nodups )

write.xlsx(list_of_datasets, 'Growth_resistance_summaryStats_NODUPS.xlsx', colNames = T, rowNames = T) 






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

# it seems that it has way less outliers than the first version, 



# GROWTH RATES ------------------------------------------------------------

# what if the max growth rate would be a better proxy for bacterial growth rather than bac growth per se
# let's try this by reading 


rates = read_csv('Summary.csv', quote = "\"") %>%
  rename(rate = `595nm_dt_Max`) %>% # growth rates
  mutate(Strain = as.character(Strain),
         Row = str_match_all(Well,'[:digit:]{1,}'), #Get plate row names from well name. 
         Col = str_match_all(Well,'[:alpha:]{1,}'), #Get plate column names from well name
         Row = factor(Row, levels = 1:12), #Make them categorical variable with set order them
         Col = factor(Col, levels = LETTERS[1:8]),
         Strain = as.factor(Strain)) %>% #Change Type column coding
  unite(ID, Strain, Plate, Well, remove = F) %>%
  select(Strain, ID, Replicate, Metformin_mM, Replicate, Plate, Well, Row, Col, rate)


# filter empty wells (some of them have growth)
rates = rates %>% filter(Strain != 'EMPTY')


# # # Fixing PG6 for rates

# according to Jen, 50 and 0 mM in PG6 are switched
# this code fixes it
sub50 = rates %>% 
  filter(Plate == 6 & Metformin_mM == 0 & Replicate == 1) %>% 
  mutate(Metformin_mM = 50)
sub50 = rates %>% 
  filter(Plate == 6 & Metformin_mM == 50 & Replicate == 1) %>% 
  mutate(Metformin_mM = 0) %>% 
  bind_rows(sub50)

rates = rates %>% 
  filter(!(Plate == 6 & Replicate == 1 & Metformin_mM %in% c(0,50)))  %>% 
  bind_rows(sub50)

# filter out the evo strains
rates = rates %>% filter(!(Strain %in% exp_str))

### duplicates ###
dups = rates %>%
  group_by(Strain) %>%
  summarise(N = n()/12) %>%
  arrange(desc(N)) %>%
  drop_na(Strain) %>%
  filter(N > 1) %>% select(Strain) %>% t %>% as.character

num = 1
rates %>% 
  filter(Strain == dups[num]) %>% 
  mutate(Metformin_mM = as.factor(Metformin_mM)) %>% 
  ggplot(aes(x = Metformin_mM, y = rate, group = Metformin_mM, fill = Metformin_mM)) + 
  geom_boxplot() +
  facet_wrap(~Plate)

rates %>% 
  filter(Strain == dups[num]) %>% 
  mutate(Metformin_mM = as.factor(Metformin_mM)) %>% 
  ggplot(aes(x = Metformin_mM, y = rate, group = Metformin_mM, fill = Metformin_mM)) + 
  geom_boxplot() 


rates


# rates summary -----------------------------------------------------------

rates.sum = rates %>% 
  filter(!Strain %in% exp_str) %>% 
  filter(!(Strain == 'NT12128' & Replicate == 2 & Metformin_mM == 50)) %>%  # this rep is bad
  filter(!Strain == 'NT12335') %>% 
  group_by(Strain, Metformin_mM) %>% 
  summarise(Mean = mean(rate, na.rm = TRUE),
            SD = sd(rate, na.rm = TRUE),
            Median = median(rate, na.rm = TRUE),
            SEM = SD/sqrt(n()))

# check how it looks 
rates.sum %>% 
  filter(Metformin_mM == 100) %>% 
  ungroup %>% 
  arrange(desc(Mean)) %>% 
  mutate(Strain = factor(Strain, levels = Strain)) %>% 
  ggplot(aes(x = Strain, y = Mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM))

# data transformation
rates.sum = rates.sum %>%
  select(-c(SD, Median, SEM)) %>% 
  pivot_wider(names_from = Metformin_mM, values_from = Mean, names_prefix = 'Metf_') %>%
  mutate(Total = Metf_0 + Metf_50 + Metf_100 + Metf_200,
         AUC50 = Metf_50/Metf_0,
         AUC100 = Metf_100/Metf_0,
         AUC200 = Metf_200/Metf_0,
         Mean = AUC50 + AUC100 + AUC200) %>% 
  select(Strain,Total, AUC50,AUC100,AUC200, Mean) 


# save info
list_of_datasets = list('unweighted' = rates.sum )

write.xlsx(list_of_datasets, 'Growth_rates_NODUPS.xlsx', colNames = T, rowNames = T) 

