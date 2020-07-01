
# libraries
library(tidyverse)
library(readxl)
library(broom)
library(openxlsx)

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


straindb = read_xlsx('strain_db.xlsx', sheet = 'strain_db')

exp_str = straindb %>% filter(Broadphenotype == 'Evolutionexperiment') %>% select(ID) %>% t %>% as.vector


# cummulative AUCs --------------------------------------------------------


bad = data %>% 
  filter(AUC_raw < 0.05, Metformin_mM == 0) %>%  arrange(ID)

w1 = 1
w2 = 
# long transformation of data
data.sum = data %>% 
  anti_join(bad) %>% 
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
  filter(Mean < 15) %>% 
  group_by(Measure) %>% 
  arrange(Mean) %>%
  mutate(ID = factor(ID, levels = ID)) %>% 
  ggplot(aes(x = ID, y = Mean)) +
  # geom_point() +
  geom_pointrange(aes(ymin = Mean - SD, ymax = Mean + SD)) +
  facet_wrap(~Measure, ncol = 1)



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



# save PCA info
list_of_datasets = list('unweighted' = data.sum %>% filter(Measure == 'uncsum'), 
                        'weighted' = data.sum %>% filter(Measure == 'wcsum'))

write.xlsx(list_of_datasets, 'Growth_resistance_summaryStats.xlsx', colNames = T, rowNames = T) 




