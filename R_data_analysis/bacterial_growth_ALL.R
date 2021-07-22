# script to analyse, and merge, datasets from bacterial growth coming from AUS and ECOREF datasets


# libraries ---------------------------------------------------------------

library(tidyverse)
library(readxl)
library(broom)
library(openxlsx)
library(here)

theme_set(theme_classic() +
            theme(axis.text=element_text(size=15),
                  axis.title=element_text(size=18,face="bold")))



# load datasets -----------------------------------------------------------

## Here, I'll load the datasets directly from their source. This way, if there
# are any modifications to those, I can always run again these lines and get 
# the latest version of them

# ECOREF data

ecoref = read_excel("D:/MRC_Postdoc/Pangenomic/pangenome_analysis/ECOREF/bacterial_growth_ECOREF/Growth_resistance_summaryStats_NODUPS_ECOREF.xlsx", 
         sheet = "unweighted") %>% 
  select(-`...1`, -AUC50, -AUC100, -AUC200,  -ID, Bact_score_mean = Mean)

# AUS data

aus = read_excel("D:/MRC_Postdoc/Pangenomic/pangenome_analysis/AUS/bacterial_growth/Growth_resistance_summaryStats_AUS.xlsx", 
      sheet = "unweighted") %>% 
  select(-`...1`, -Measure, -Bact_score_sd)

# metadata
metadata = read_excel("D:/MRC_Postdoc/Pangenomic/pangenome_analysis/metadata/MAIN_metadata.xlsx", 
                      sheet = "metadata")



### join both datasets

# first check that they share the same column names


names(ecoref) %in% names(aus)
names(aus) %in% names(ecoref)

data = ecoref %>% select(names(aus)) %>% 
  bind_rows(aus) %>% 
  select(-Plate) %>% 
  # filter(phylogroup != 'Non Escherichia') %>% 
  left_join(metadata %>% select(Strain = ID, everything()))


write.csv(data, 'bacterial_growth_ALL.csv')

data %>% 
  filter(phylogroup != 'Non Escherichia') %>%
  filter(!(Strain %in% c('NT12335','NT12332'))) %>% 
  ggplot(aes(x = fct_reorder(Strain, Bact_score_mean), y = Bact_score_mean, 
             color = Origin)) +
  geom_point(size = 2, alpha = 0.5)

ggsave(here('summary', 'bact_scores.pdf'))


data %>% 
  filter(phylogroup != 'Non Escherichia') %>%
  filter(!(Strain %in% c('NT12335','NT12332'))) %>% 
  ggplot(aes(x = Origin, y = Bact_score_mean, 
             fill = Origin)) +
  geom_boxplot()

ggsave(here('summary', 'boxplot_bact_scores.pdf'))



