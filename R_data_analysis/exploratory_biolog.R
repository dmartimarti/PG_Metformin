
# libraries
library(tidyverse)
library(readxl)
library(broom)

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

data %>% filter(Strain != 'EMPTY')

# data = read.table('data_rawAUC_PG.txt', sep = '\t', header = TRUE)
#data = as_tibble(data) %>% 
#  unite(ID, Strain, Plate, Well, remove = F) 


met0 = data %>% filter(Metformin_mM == 0) %>%
  select(ID, Replicate, AUC_raw) %>%
  pivot_wider(names_from = Replicate, values_from = AUC_raw, names_prefix = 'rep_')

weird_st = met0 %>% filter(rep_1 < 0.05 | rep_2 < 0.05 | rep_3 < 0.05) %>%
  select(ID) %>% as.vector %>% t 


met0 %>% filter(!ID %in% weird_st) %>%
  ggplot(aes(x = rep_1, y = rep_2)) +
  geom_point() + 
  theme_classic()


met0  %>%
  ggplot(aes(x = rep_1, y = rep_2)) +
  geom_point() + 
  theme_classic()


met0 %>% ggplot(aes(x = rep_1, y = rep_3)) +
  geom_point() + 
  theme_classic()

met0 %>% ggplot(aes(x = rep_2, y = rep_3)) +
  geom_point() + 
  theme_classic()

model = lm(rep_1 ~ rep_2, data = met0)
summary(model)
model = lm(rep_1 ~ rep_3, data = met0)
summary(model)
model = lm(rep_2 ~ rep_3, data = met0)
summary(model)


# this is telling us something weird is happening
# lets remove evolution strains

straindb = read_xlsx('strain_db.xlsx', sheet = 'strain_db')

exp_str = straindb %>% filter(Broadphenotype == 'Evolutionexperiment') %>% select(ID) %>% t %>% as.vector

met0_filt = data %>% filter(Metformin_mM == 0) %>%
  filter(!Strain %in% exp_str) %>%
  select(ID, Replicate, AUC_raw) %>%
  pivot_wider(names_from = Replicate, values_from = AUC_raw, names_prefix = 'rep_')

# remove weird strains with low growth
# weird_st = met0_filt %>% filter(rep_1 < 0.05 | rep_2 < 0.05 | rep_3 < 0.05) %>%
#   select(ID) %>% as.vector %>% t 

# met0_filt = met0_filt %>% filter(!ID %in% weird_st)

met0_filt  %>%
  ggplot(aes(x = rep_1, y = rep_2)) +
  geom_point() + 
  theme_classic()


met0_filt %>% ggplot(aes(x = rep_1, y = rep_3)) +
  geom_point() + 
  theme_classic()

met0_filt %>% ggplot(aes(x = rep_2, y = rep_3)) +
  geom_point() + 
  theme_classic()

model = glance(lm(rep_1 ~ rep_2, data = met0_filt))
rep1_2 = glance(summary(model))
model = lm(rep_1 ~ rep_3, data = met0_filt)
rep1_3 = glance(summary(model))
model = lm(rep_2 ~ rep_3, data = met0_filt)
rep2_3 = glance(summary(model))

stats = data.frame()
stats = rbind(rep1_2, rep1_3, rep2_3)




# calculate pairwise lm functions
rep.lm = function(df, met = 0, thr = 0){
  # filter by values
  weird_st = df %>% filter(rep_1 <= thr | rep_2 <= thr | rep_3 <= thr) %>%
    select(ID) %>% as.vector %>% t 
  df = df %>% filter(!ID %in% weird_st)
  
  rep1_2 = glance(summary(lm(rep_1 ~ rep_2, data = df)))
  rep1_3 = glance(summary(lm(rep_1 ~ rep_3, data = df)))
  rep2_3 = glance(summary(lm(rep_2 ~ rep_3, data = df)))
  stats = data.frame()
  stats = rbind(rep1_2, rep1_3, rep2_3)
  stats['Metformin_mM'] = met
  stats['Comparison'] = c('rep1_2', 'rep1_3', 'rep2_3')
  return(stats %>% select(Comparison, Metformin_mM, everything()))
}

met0_filt = data %>% filter(Metformin_mM == 0) %>%
  filter(!Strain %in% exp_str) %>%
  select(ID, Replicate, AUC_raw) %>%
  pivot_wider(names_from = Replicate, values_from = AUC_raw, names_prefix = 'rep_')
met50_filt = data %>% filter(Metformin_mM == 50) %>%
  filter(!Strain %in% exp_str) %>%
  select(ID, Replicate, AUC_raw) %>%
  pivot_wider(names_from = Replicate, values_from = AUC_raw, names_prefix = 'rep_')
met100_filt = data %>% filter(Metformin_mM == 100) %>%
  filter(!Strain %in% exp_str) %>%
  select(ID, Replicate, AUC_raw) %>%
  pivot_wider(names_from = Replicate, values_from = AUC_raw, names_prefix = 'rep_')
met200_filt = data %>% filter(Metformin_mM == 200) %>%
  filter(!Strain %in% exp_str) %>%
  select(ID, Replicate, AUC_raw) %>%
  pivot_wider(names_from = Replicate, values_from = AUC_raw, names_prefix = 'rep_')


met0_stats    =  rep.lm(met0_filt, thr = 0.1, met = 0)
met50_stats =   rep.lm(met50_filt,thr = 0.1, met = 50)
met100_stats = rep.lm(met100_filt, thr = 0.1,met = 100)
met200_stats = rep.lm(met200_filt,thr = 0.1, met = 200)

data_fits = rbind(met0_stats, met50_stats, met100_stats, met200_stats)


data_fits %>%
  mutate(Metformin_mM = as.factor(Metformin_mM)) %>%
  ggplot(aes(x = Comparison, y = adj.r.squared, colour = Metformin_mM)) +
  geom_point(size = 4) +
  theme_classic()




