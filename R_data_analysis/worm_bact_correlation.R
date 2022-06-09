
# libraries ---------------------------------------------------------------

library(tidyverse)
library(readxl)
library(readr)
library(here)
library(cowplot)
library(ggpubr)


# read the data -----------------------------------------------------------

metadata = read_excel("~/Documents/MRC_postdoc/Pangenomic/metadata/MAIN_metadata.xlsx")

worm_pheno = 
  read_excel("~/Documents/MRC_postdoc/Pangenomic/pangenome_analysis/ALL/worm_imaging/ALL_worm_FC.xlsx")

bac_pheno = 
  read_csv("~/Documents/MRC_postdoc/Pangenomic/pangenome_analysis/ALL/bacterial_growth/bacterial_growth_ALL.csv")


# filter and merge datasets -----------------------------------------------

discard = c('cladeI','Non Escherichia','E or cladeI',
            'cladeII','cladeIII','cladeIV','cladeV',
            'albertii','fergusonii')

metadata = metadata %>% 
  filter(Discard == 'No',
         !(phylogroup %in% discard),
         Origin %in% c('AUS','ECOREF')) %>% 
  distinct(.keep_all = T) %>% 
  drop_na(Broadphenotype)

# filter worm phenotype
worm_pheno = worm_pheno %>% filter(ID %in% metadata$ID)

# filter bact phenotype
bac_pheno = bac_pheno %>% filter(Strain %in% metadata$ID)


# join both phenotypes
phenotypes = bac_pheno %>% 
  select(ID = Strain, bact_pheno = Bact_score_mean) %>% 
  left_join(worm_pheno %>% 
              select(ID, worm_pheno = Mean_FC)) %>% 
  left_join(metadata %>% 
              select(ID, Broadphenotype, phylogroup))


# plot correlation --------------------------------------------------------

phenotypes %>% 
  filter(bact_pheno < 4) %>% 
  ggplot(aes(x = bact_pheno, y = worm_pheno)) +
  geom_smooth(method = 'lm') +
  geom_point() +
  theme_cowplot(15)


phenotypes %>% 
  filter(bact_pheno < 4) %>% 
  filter(Broadphenotype != 'Unknown') %>% 
  ggplot(aes(x = bact_pheno, y = worm_pheno)) +
  geom_smooth(method = 'lm') +
  geom_point(alpha = 0.5) +
  stat_cor(method = "pearson", label.x = 0.6, label.y = 0.1) +
  labs(
    y = 'Worm phenotype\n(Mean fold change)',
    x = 'Bacterial phenotype\n(Mean resistance score)'
  ) +
  scale_x_continuous(breaks = seq(0.5, 2.5, 1)) +
  facet_wrap(~Broadphenotype) +
  theme_cowplot(15) 


ggsave('exploration/correlation_phenotype.pdf',
       height = 5, width = 9)



phenotypes %>% 
  filter(bact_pheno < 4) %>% 
  filter(Broadphenotype != 'Unknown') %>% 
  ggplot(aes(x = bact_pheno, y = worm_pheno)) +
  geom_smooth(method = 'lm') +
  geom_point() +
  stat_cor(method = "pearson", label.x = 0.6, label.y = 0.1) +
  labs(
    y = 'Worm phenotype\n(Mean fold change)',
    x = 'Bacterial phenotype\n(Mean resistance score)'
  ) +
  facet_wrap(~phylogroup) +
  theme_cowplot(15) 



