
# libraries ---------------------------------------------------------------

library(tidyverse)
library(openxlsx)
library(readxl)



# datasets ----------------------------------------------------------------

# AUS info
AUS_biofilm = read_excel("201201_biofilm_annotation_australian_strains_AVERAGE.xlsx")

# ECOREF biofilm
ecoref_bio = read_excel("20909_biofilm_strain_annotation_ECOREF.xlsx", 
                                                      sheet = "rep4")

# Main metadata
metadata = read_excel("MAIN_metadata.xlsx", 
                            sheet = "metadata")
# phylogroup analysis
phylog = read_excel("analysis_2021-04-29_191026_phylogroups.xlsx") %>% 
  select(-`...2`,-`...3`,-`...4`,-`...8`,-`...9`) %>% 
  rename(cont = `...7`)

# quast_seqs
quast_seqs = read_excel("assembly_report_quast.xlsx", 
                                    sheet = "transposed_report")



# the objective here is to merge the new entries into the MAIN metadata


# operations --------------------------------------------------------------

### ECOREF ####

PG_missing = metadata %>% 
  filter(is.na(Assembly)) 

# filter ECOREF, modify its assembly name
quast_ECOREF = quast_seqs %>% 
  filter(Origin == 'ECOREF') %>% 
  select(Assembly, strain, Origin, REPEAT) %>% 
  filter(strain %in% PG_missing$ID) %>%  # filter repeated strains
  mutate(Assembly = paste0(strain,'.fasta')) 

# filter phylogroup
phylog_ECOREF = phylog %>% 
  filter(Genome %in% quast_ECOREF$Assembly) %>% 
  mutate(phylogroup = case_when(str_detect(phylogroup, 'Unknown') ~ mash_group,
                                TRUE ~ phylogroup)) %>% 
  rename(Assembly = Genome) %>% 
  select(Assembly, phylogroup, mash_group)

# join both datasets
missing_ECOREF = phylog_ECOREF %>% left_join(quast_ECOREF) %>% 
  select(-REPEAT) %>% 
  rename(ID = strain,
         fasta = Assembly) %>% 
  mutate(Assembly = ID)

# remove columns from original metadata and substitute them for new info

PG_missing %>%
  select(-Assembly, phylogroup, mash_group) %>% 
  left_join(missing_ECOREF)


# save list
new_PG = PG_missing %>% 
  filter((ID %in% missing_ECOREF$ID)) %>% 
  select(-Assembly, -phylogroup, -mash_group, -fasta) %>% 
  left_join(missing_ECOREF) %>% 
  select(names(PG_missing))


# PG_missing %>% 
#   filter(!(ID %in% missing_ECOREF$ID)) %>% view
# there are two strains that didn't get sequenced this time but we don't need them
# one is a ev. experiment strain, and another is a strain with a 
# similar assembly as one of the strains contained here



# remove old strains and add new

metadata = metadata %>% 
  filter(!(ID %in% new_PG$ID)) %>% 
  mutate(Assembly = as.character(Assembly)) %>% 
  bind_rows(new_PG) %>% 
  mutate(Origin = 'ECOREF')


### AUS ####

AUS_biofilm = AUS_biofilm %>% 
  mutate(PG = paste0('AUS_',PG)) %>% 
  mutate(ID = Strain,
         Genome = paste0(Strain,'.fasta')) %>% 
  select(-Strain)  


phylog %>% 
  filter(Genome %in% AUS_biofilm$Genome)












### remove contaminations ####






