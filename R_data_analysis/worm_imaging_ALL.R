# script to analyse, and merge, datasets from worm acs-2 assays 
# coming from AUS and ECOREF datasets


# libraries ---------------------------------------------------------------

library(tidyverse)
library(readr)
library(readxl)
library(broom)
library(openxlsx)
library(here)

theme_set(theme_classic() +
            theme(axis.text=element_text(size=15),
                  axis.title=element_text(size=18,face="bold")))



# read datasets -----------------------------------------------------------

# read files containing the FC of worm brightness

# AUS strains
AUS_FC = read_csv("D:/MRC_Postdoc/Pangenomic/pangenome_analysis/AUS/worm_imaging/analysis/FC_means_unique.csv") %>% 
  select(-PG) %>% 
  mutate(ID = as.factor(ID))

ECO_FC = read_csv("D:/MRC_Postdoc/Pangenomic/pangenome_analysis/ECOREF/worm_imaging_ECOREF/analysis/FC_means_unique.csv") %>% 
  select(-Strainname, -PG) %>% 
  mutate(ID = as.factor(ID))


metadata = read_excel("D:/MRC_Postdoc/Pangenomic/pangenome_analysis/metadata/MAIN_metadata.xlsx", 
                            sheet = "metadata")


AUS_FC

ECO_FC

all_FC = ECO_FC %>% bind_rows(AUS_FC)

# check that we really have only one value per strain
length(unique(all_FC$ID)) == dim(all_FC)[1]

all_FC_metadata = all_FC %>% left_join(metadata)

all_FC_metadata %>% 
  drop_na(Mean_FC) %>% 
  ggplot(aes(x = fct_reorder(ID, Mean_FC), y = Mean_FC, color = Origin)) +
  geom_point(size = 2, alpha = 0.7)


all_FC_metadata %>% 
  drop_na(Mean_FC) %>% 
  filter(Annotation_50mM == 'normal') %>% 
  ggplot(aes(x = fct_reorder(ID, Mean_FC), y = Mean_FC, color = Origin)) +
  geom_point(size = 2, alpha = 0.7)


all_FC_metadata %>% 
  drop_na(Mean_FC) %>% 
  filter(Annotation_50mM == 'normal') %>% 
  mutate(fasta = str_sub(fasta,1, -7)) %>% 
  drop_na(fasta) %>% 
  select(IDs = fasta, FC_worm = Mean_FC) %>% 
  write_delim('worm_phenotype_no_biofilm.txt', delim = '\t')





