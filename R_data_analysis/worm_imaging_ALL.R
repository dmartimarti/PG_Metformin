# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# script to analyse, and merge, datasets from worm acs-2 assays 
# coming from AUS and ECOREF datasets
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# libraries ---------------------------------------------------------------

library(tidyverse)
library(readr)
library(readxl)
library(broom)
library(openxlsx)
library(here)
library(cowplot)

theme_set(theme_classic() +
            theme(axis.text=element_text(size=15),
                  axis.title=element_text(size=18,face="bold")))



# read datasets -----------------------------------------------------------

# read files containing the FC of worm brightness

# AUS strains
AUS_FC = read_csv("D:/MRC_Postdoc/Pangenomic/pangenome_analysis/AUS/worm_imaging/analysis/FC_means_unique.csv") %>% 
  select(-PG) %>% 
  mutate(ID = as.factor(ID))

# ECOREF strains
ECO_FC = read_csv("D:/MRC_Postdoc/Pangenomic/pangenome_analysis/ECOREF/worm_imaging_ECOREF/analysis/FC_means_unique.csv") %>% 
  select(-Strainname, -PG) %>% 
  mutate(ID = as.factor(ID))


metadata = read_excel("D:/MRC_Postdoc/Pangenomic/pangenome_analysis/metadata/MAIN_metadata.xlsx", 
                            sheet = "metadata")

# for Mac
metadata = read_excel(here("../../../metadata/","MAIN_metadata.xlsx"), 
                      sheet = "metadata")

AUS_FC

ECO_FC

all_FC = ECO_FC %>% bind_rows(AUS_FC)

# check that we really have only one value per strain
length(unique(all_FC$ID)) == dim(all_FC)[1]

all_FC_metadata = all_FC %>% 
  select(-Annotation_0mM, -Annotation_50mM) %>% 
  left_join(metadata, by = c('ID'))

all_FC_metadata %>% 
  drop_na(Mean_FC) %>% 
  ggplot(aes(x = fct_reorder(ID, desc(Mean_FC)), y = Mean_FC, color = Origin)) +
  geom_point(size = 2, alpha = 0.7) +
  theme(axis.text.x = element_blank())


all_FC_metadata %>% 
  drop_na(Mean_FC) %>% 
  distinct(ID, .keep_all = T) %>% 
  ggplot(aes(x = fct_reorder(ID, desc(Mean_FC)), y = Mean_FC)) +
  geom_point(size = 1, alpha = 0.7) +
  geom_errorbar(aes(ymax = Mean_FC + SD_FC, ymin = Mean_FC - SD_FC), size = 0.2) +
  labs(
    x = 'Strains',
    y = 'Mean FC (± SD)'
  ) + 
  theme_cowplot(17) +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank())

ggsave(here('exploration', 'worm_FC_overview.pdf'),
       height = 7, width = 9)

all_FC_metadata %>% 
  drop_na(Mean_FC) %>% 
  filter(Annotation_50mM == 'normal') %>% 
  ggplot(aes(x = fct_reorder(ID, Mean_FC), y = Mean_FC, color = Origin)) +
  geom_point(size = 2, alpha = 0.7)


#### save merged data ####

all_FC_metadata %>% 
  write.xlsx('ALL_worm_FC.xlsx', sheetName = 'FC_per_strain',
             overwrite = T)



removals = c('Non Escherichia', 'cladeI', 'cladeII', 'cladeIII', 
             'cladeIV', 'cladeV', 'albertii', 'fergusonii')



# datasets with stats -----------------------------------------------------

# AUS strains
AUS_stats = read_excel("D:/MRC_Postdoc/Pangenomic/pangenome_analysis/AUS/worm_imaging/analysis/worm_imaging_stats.xlsx",
                      sheet = 'Stats_per_plate') %>% 
  select(-ID) %>% 
  left_join(read_csv("D:/MRC_Postdoc/Pangenomic/pangenome_analysis/AUS/worm_imaging/analysis/FC_means_unique.csv") %>% 
              select(PG:Annotation_50mM)) %>% 
  select(-PG) %>% 
  mutate(ID = as.factor(ID))
  
# ECOREF strains
ECO_stats = read_excel("D:/MRC_Postdoc/Pangenomic/pangenome_analysis/ECOREF/worm_imaging_ECOREF/analysis/worm_imaging_stats.xlsx",
                     sheet = 'Stats_per_plate')   %>% 
  select(-Strainname, -PG) %>% 
  mutate(ID = as.factor(ID))
  
AUS_stats

ECO_stats

all_stats = ECO_stats %>% 
  select(-phylogroup) %>% 
  bind_rows(AUS_stats) %>% 
  drop_na(ID) %>% 
  distinct(ID, .keep_all = T)

# check that we really have only one value per strain
length(unique(all_stats$ID)) == dim(all_stats)[1]


# join datasets into one that has all the relevant info
all_stats_metadata = all_stats %>% 
  select(Well, ID, FC:FDR_stars) %>% 
  left_join(metadata, by = c('ID', 'Well')) 


#### save stats merged ####
all_stats_metadata %>% 
  write.xlsx('ALL_worm_stats.xlsx', sheetName = 'Stats_per_strain',
             overwrite = T)




# worm datasets for Pyseer -----------------------------------------------------



# save all genomes without biofilm at 50mM
all_FC_metadata %>% 
  drop_na(Mean_FC) %>% 
  filter(Annotation_50mM == 'normal') %>% 
  mutate(fasta = str_sub(fasta,1, -7)) %>% 
  drop_na(fasta) %>% 
  select(IDs = fasta, FC_worm = Mean_FC) %>% 
  write_delim('worm_phenotype_no_biofilm.txt', delim = '\t')


# save all genomes 
all_FC_metadata %>% 
  drop_na(Mean_FC) %>% 
  # filter(Annotation_50mM == 'normal') %>% 
  mutate(fasta = str_sub(fasta,1, -7)) %>% 
  drop_na(fasta) %>% 
  select(IDs = fasta, FC_worm = Mean_FC) %>% 
  write_delim('worm_phenotype_ALL.txt', delim = '\t')

# save ECO genomes without biofilm
all_FC_metadata %>% 
  drop_na(Mean_FC) %>% 
  filter(Origin == 'ECOREF') %>% 
  filter(Annotation_50mM == 'normal') %>% 
  mutate(fasta = str_sub(fasta,1, -7)) %>% 
  drop_na(fasta) %>% 
  select(IDs = fasta, FC_worm = Mean_FC) %>% 
  write_delim('worm_phenotype_ECOREF_no_biofilm.txt', delim = '\t')

# save ECO genomes
all_FC_metadata %>% 
  drop_na(Mean_FC) %>% 
  filter(Origin == 'ECOREF') %>% 
  # filter(Annotation_50mM == 'normal') %>% 
  mutate(fasta = str_sub(fasta,1, -7)) %>% 
  drop_na(fasta) %>% 
  select(IDs = fasta, FC_worm = Mean_FC) %>% 
  write_delim('worm_phenotype_ECOREF.txt', delim = '\t')


# save AUS genomes without biofilm
all_FC_metadata %>% 
  drop_na(Mean_FC) %>% 
  filter(Origin == 'AUS') %>% 
  filter(Annotation_50mM == 'normal') %>% 
  mutate(fasta = str_sub(fasta,1, -7)) %>% 
  drop_na(fasta) %>% 
  select(IDs = fasta, FC_worm = Mean_FC) %>% 
  write_delim('worm_phenotype_AUS_no_biofilm.txt', delim = '\t')

# save AUS genomes
all_FC_metadata %>% 
  drop_na(Mean_FC) %>% 
  filter(Origin == 'AUS') %>% 
  # filter(Annotation_50mM == 'normal') %>% 
  mutate(fasta = str_sub(fasta,1, -7)) %>% 
  drop_na(fasta) %>% 
  select(IDs = fasta, FC_worm = Mean_FC) %>% 
  write_delim('worm_phenotype_AUS.txt', delim = '\t')






# biofilm producers (for Jen) ---------------------------------------------

# Jen asked me to pull the strains that produce superbio and do pyseer to them

library(readxl)
biofilm = read_excel("metf_induced_biofilm.xlsx")


biofilm_strains = biofilm %>% pull(Strain)

# save all genomes 
### THIS VERSION IS WRONG
# I took the strains that were producing biofilm and analysed them, but
# they should be labelled as biofilm (1) or non-biofilm (2)
# all_FC_metadata %>% 
#   filter(ID %in% biofilm_strains) %>% 
#   drop_na(Mean_FC) %>% 
#   # filter(Annotation_50mM == 'normal') %>% 
#   mutate(fasta = str_sub(fasta,1, -7)) %>% 
#   mutate(fasta = case_when(is.na(fasta) ~ ID,
#                            TRUE ~ fasta)) %>% 
#   drop_na(fasta) %>% 
#   select(IDs = fasta, FC_worm = Mean_FC) %>% 
#   write_delim('worm_phenotype_normal2superbio.txt', delim = '\t')


all_FC_metadata %>%
  # filter(ID %in% biofilm_strains) %>%
  drop_na(Mean_FC) %>%
  # filter(Annotation_50mM == 'normal') %>%
  mutate(fasta = str_sub(fasta,1, -7)) %>%
  mutate(fasta = case_when(is.na(fasta) ~ ID,
                           TRUE ~ fasta)) %>%
  drop_na(fasta) %>%
  select(IDs = fasta, FC_worm = Mean_FC) %>%
  mutate(Biofilm = case_when(IDs %in% biofilm_strains ~ 1,
                             TRUE ~ 2)) %>% 
  write_delim('worm_phenotype_normal2superbio.txt', delim = '\t')



# this chunk gets the bacteria phenotype as by Jen's biofilm data, 
# but only taking the biofilm info (not bacterial growth or worm phenotype)
all_FC_metadata %>%
  filter(Annotation_0mM == 'normal') %>%
  mutate(fasta = str_sub(fasta,1, -7)) %>%
  mutate(fasta = case_when(is.na(fasta) ~ ID,
                           TRUE ~ fasta)) %>%
  drop_na(fasta) %>%
  select(IDs = fasta) %>%
  mutate(Biofilm = case_when(IDs %in% biofilm_strains ~ 1,
                             TRUE ~ 2)) %>% 
  write_delim('PG_bacterial_normal2superbio.txt', delim = '\t')


# save all genomes 
all_FC_metadata %>% 
  filter(ID %in% biofilm_strains) %>% 
  drop_na(Mean_FC) %>% 
  # filter(Annotation_50mM == 'normal') %>% 
  mutate(fasta = str_sub(fasta,1, -7)) %>% 
  # drop_na(fasta) %>% 
  # select(IDs = fasta, FC_worm = Mean_FC)  %>% 
  view





all_FC_metadata %>%
  mutate(Path = paste0('/rds/general/user/dmarti14/home/pangenome_study/complete/assemblies/no_evo/',fasta)) %>% 
  filter(Annotation_0mM == 'normal') %>%
  filter(!(phylogroup %in% removals)) %>% 
  filter(Discard == 'No') %>% 
  mutate(fasta = str_sub(fasta,1, -7)) %>%
  mutate(fasta = case_when(is.na(fasta) ~ ID,
                           TRUE ~ fasta)) %>%
  drop_na(fasta) %>%
  select(ID = fasta, Path) %>%
  mutate(Phenotype = case_when(ID %in% biofilm_strains ~ 1,
                             TRUE ~ 0)) %>% 
  select(ID, Phenotype, Path)  %>% 
  distinct(ID, .keep_all = TRUE) %>% 
  write_delim('dbgwas_biofilm_PG.txt', delim = '\t')
  


# DBGWAS versions ---------------------------------------------------------

# save all genomes without biofilm at 50mM
all_FC_metadata %>% 
  mutate(Path = paste0('/rds/general/user/dmarti14/home/pangenome_study/complete/assemblies/no_evo/',fasta)) %>% 
  drop_na(Mean_FC) %>% 
  filter(Annotation_50mM == 'normal') %>% 
  filter(!(phylogroup %in% removals)) %>% 
  filter(Discard == 'No') %>% 
  mutate(fasta = str_sub(fasta,1, -7)) %>% 
  drop_na(fasta) %>% 
  select(ID = fasta, Phenotype = Mean_FC, Path) %>% 
  distinct(ID, .keep_all = TRUE) %>% 
  write_delim('worm_phenotype_no_biofilm_DBGWAS.txt', delim = '\t')


# save all genomes 
all_FC_metadata %>% 
  mutate(Path = paste0('/rds/general/user/dmarti14/home/pangenome_study/complete/assemblies/no_evo/',fasta)) %>% 
  drop_na(Mean_FC) %>% 
  filter(!(phylogroup %in% removals)) %>% 
  filter(Discard == 'No') %>% 
  mutate(fasta = str_sub(fasta,1, -7)) %>% 
  drop_na(fasta) %>% 
  select(ID = fasta, Phenotype = Mean_FC, Path) %>% 
  distinct(ID, .keep_all = TRUE) %>% 
  write_delim('worm_phenotype_ALL_DBGWAS.txt', delim = '\t')






