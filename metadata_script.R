
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
metadata = read_excel("MAIN_metadata_OLD.xlsx", 
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
  select(Assembly, strain, Origin, REPEAT, Notes = `Dani's notes`) %>% 
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
  select(-Assembly, -phylogroup,-Notes, -mash_group, -fasta) %>% 
  left_join(missing_ECOREF) %>% 
  select(names(PG_missing))


# PG_missing %>% 
#   filter(!(ID %in% missing_ECOREF$ID)) %>% view
# there are two strains that didn't get sequenced this time but we don't need them
# one is a ev. experiment strain, and another is a strain with a 
# similar assembly as one of the strains contained here



# remove old strains and add new

metadata_ECOREF = metadata %>% 
  filter(!(ID %in% new_PG$ID)) %>% 
  mutate(Assembly = as.character(Assembly)) %>% 
  bind_rows(new_PG) %>% 
  mutate(Origin = 'ECOREF')


# add data about biofilm formation
metadata_ECOREF = metadata_ECOREF %>% 
  left_join(ecoref_bio %>% 
              rename(Annotation_0mM = `0mM_annotation`,
                     Annotation_50mM = `50mM_annotation`) %>% 
              select(ID = Strain, Annotation_0mM,Annotation_50mM, Well)) %>% 
  select(-F_annotation, -J_annotation)



### AUS ####

AUS_biofilm = AUS_biofilm %>% 
  mutate(PG = paste0('AUS_',PG)) %>% 
  mutate(ID = Strain,
         Genome = paste0(Strain,'.fasta')) %>% 
  select(-Strain)  


AUS_missing = phylog %>% 
  filter(Genome %in% AUS_biofilm$Genome) %>% 
  left_join(AUS_biofilm) %>% 
  rename(fasta = Genome) %>% 
  select(-Position) %>% 
  mutate(ID = as.character(ID),
         Origin = 'AUS')



## JOIN BOTH DATASETS
metadata_ECOREF_AUS = metadata_ECOREF %>% 
  bind_rows(AUS_missing)


### Others ####

quast_others = quast_seqs %>% 
  filter(Origin %in% c('OP50', 'Germany', 'ECOREF_Ev')) %>% 
  select(strain, Origin, Notes = `Dani's notes`) %>% 
  mutate(strain = case_when(strain == '2.2999999999999998' ~ '2.3',
                            TRUE ~ strain),
         Genome = paste0(strain, '.fasta'))
# join and change
others = phylog %>% 
  filter(Genome %in% quast_others$Genome) %>% 
  select(Genome, phylogroup) %>% 
  left_join(quast_others) %>% 
  rename(fasta = Genome,
         ID = strain)


### JOIN ALL DATASETS INTO ONE ####

full_metadata = metadata_ECOREF_AUS %>% 
  bind_rows(others) %>% 
  rename(AUS_notes = notes)


unique(full_metadata$Annotation_50mM)
unique(full_metadata$Annotation_0mM)



# small tweeks
# modify the biofilm annotation

full_metadata = full_metadata %>% 
  mutate(Annotation_0mM = case_when(Annotation_0mM == 'NB' ~ 'normal',
                                    Annotation_0mM == 'super_biofilm' ~ 'super_bio',
                                    TRUE ~ Annotation_0mM),
         Annotation_50mM = case_when(Annotation_50mM == 'super_biofilm' ~ 'super_bio',
                                     TRUE ~ Annotation_50mM)) %>% 
  mutate(Plate = case_when(Origin == 'AUS' ~ str_sub(PG, start=-1),
                           Origin != 'AUS' ~ as.character(Plate))) %>% 
  mutate(Plate_Well = paste(PG,Plate,sep = '_')) 




full_metadata



list_of_datasets = list(
  'metadata' = full_metadata
)

write.xlsx(list_of_datasets, 'MAIN_metadata.xlsx')




full_metadata %>% 
  filter(!is.na(Notes))


# annotate duplicates ---------------------
library(readxl)
MAIN_metadata <- read_excel("MAIN_metadata.xlsx", 
                            sheet = "metadata")

dups = MAIN_metadata$ID[duplicated(MAIN_metadata$ID)]

full_meta_dups = MAIN_metadata %>% 
  mutate(Duplicate = case_when(ID %in% dups ~ 'Yes',
                               TRUE ~ 'No')) %>% 
  group_by(ID, Duplicate) %>% 
  mutate(Duplicate_n = seq_along(Duplicate)) %>% 
  list('metadata' = .)

write.xlsx(full_meta_dups, 'MAIN_metadata.xlsx')




# RNAseq landscape version ------------------------------------------------



MAIN_metadata = read_excel("MAIN_metadata.xlsx", 
                            sheet = "metadata")

# load the strains that produce biofilms in control conditions (LB media)
bcsA_strains =
  read_excel("RNAseq_project/E. coli collection EcoRef and Aus- bcsA strains .xlsx") %>% 
  pull(Strain)

bcsA_strains

# filter
rnaseq_meta = MAIN_metadata %>% 
  select(-Annotation_0mM, -Annotation_50mM, -Annotation) %>% 
  filter(Origin %in% c('AUS', 'ECOREF'),
         Duplicate_n == 1,
         Discard == 'No') %>% 
  filter(!(ID %in% bcsA_strains)) %>% 
  mutate(RNAseq = 'valid')


rnaseq_meta %>% 
  distinct(ID)

# save the file
write.xlsx(rnaseq_meta, 'RNAseq_project/RNAseq_metadata.xlsx')



# plate layout ------------------------------------------------------------


rnaseq_meta = read_excel("RNAseq_project/RNAseq_metadata.xlsx")

plates  = unique(rnaseq_meta$PG)

wells_expanded = expand_grid(expand_grid(col = LETTERS[1:8], row = seq(1,12)) %>% 
                               unite(Well, col, row, sep = ''), 
                             plates) %>% 
  rename(PG = plates)

plate_layout = MAIN_metadata %>% 
  left_join(rnaseq_meta) %>%
  filter(PG %in% plates) %>% 
  full_join(wells_expanded) %>% 
  mutate(NewID = case_when(RNAseq == 'valid' ~ ID,
                            TRUE ~ 'NOT VALID')) 
  



  

# initialise variables
plates = unique(rnaseq_meta$PG)
row_names = LETTERS[1:8]
col_names = seq(1,12)

plates_list = list()
for (plate in plates) {
  
  plate_name = as.character(plate)
  
  temp_ids = plate_layout %>% 
    mutate(Column = str_extract(Well, '\\w'),
           Row = str_extract(Well, '\\d{1,}')) %>% 
    filter(PG == plate) %>% 
    mutate(Row = as.numeric(Row)) %>% 
    arrange(Column, Row) %>% 
    pull(NewID)  
  
  temp_matrix = matrix(data = temp_ids, 
                       nrow = 8, ncol = 12, byrow = T)
  
  rownames(temp_matrix) = row_names
  colnames(temp_matrix) = col_names
  
  temp_list = list(plate_name = temp_matrix)
  
  plates_list = append(plates_list, temp_list)
}

# fix naming of plates in list
names(plates_list) = plates

sink("RNAseq_project/RNAseq_plate_layout.csv", type="output")
invisible(lapply(names(plates_list), 
                 function(x) { print(x)
                   dput(write.csv(plates_list[[x]])) } ))
sink()


write.xlsx(plates_list, 'RNAseq_project/RNAseq_plate_layout_multipage.xlsx',
           overwrite = T,
           rownames = T)




Metadata_phenotype_extended = read_excel("Metadata_phenotype_extended.xlsx")


Metadata_phenotype_extended %>% 
  filter(Discard == 'No',
         Origin %in% c('AUS', 'ECOREF')) %>% 
  distinct(fasta, .keep_all = T) %>% count(extended_phenotype)




Metadata_phenotype_extended %>% 
  filter(Discard == 'No',
         Origin %in% c('AUS', 'ECOREF')) %>% 
  distinct(fasta, .keep_all = T) %>% 
  filter(phylogroup != 'E or cladeI') %>% 
  count(phylogroup, extended_phenotype) %>% 
  group_by(extended_phenotype) %>% 
  mutate(prop = (n / sum(n)) * 100) %>% 
  drop_na() %>% 
  filter(extended_phenotype != 'Unknown') %>% 
  ggplot(aes(y = prop, x = phylogroup, fill = phylogroup)) +
  geom_bar(stat = 'identity', color = 'black') +
  facet_wrap(~extended_phenotype) +
  scale_fill_manual(
    values = c(
    'A' = '#EA392D',
    'B1' = '#277FEA',
    'B2' = '#0FEB3B',
    'C' = '#5910EA',
    'D' = '#EB8418',
    'E' = '#EAE21D',
    'F' = '#967944',
    'G' = '#E11CEB')
  ) +
  labs(
    y = '% of total',
    x = NULL
  ) +
  theme_cowplot(14)


