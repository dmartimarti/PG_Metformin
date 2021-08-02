library(readr)
library(tidyverse)
library(here)

theme_set(theme_classic() +
            theme(axis.text.x = element_text(size = 13, color = 'black'),
                  axis.text.y = element_text(size = 13, color = 'black'),
                  axis.title.x = element_text(face = "bold", size = 13, color = 'black'),
                  axis.title.y = element_text(face = "bold", size = 13, color = 'black')))
# Read data ---------------------------------------------------------------


muts =  read_csv("D:/MRC_Postdoc/Pangenomic/mutants_analysis/results/mutations_summary.csv")

# prep the datatable
muts = muts %>% 
  mutate(Strain = as.factor(Strain), 
         TYPE = as.factor(TYPE),
         FTYPE = as.factor(FTYPE),
         STRAND = as.factor(STRAND)) %>% 
  separate(EFFECT, into = c('Variant', 'Variant_nt','Variant_aa'), sep = ' ') %>% 
  mutate(Variant = as.factor(Variant)) %>% 
  mutate(GENE = case_when(PRODUCT == 'hypothetical protein' ~ 'hypothetical protein',
                          TRUE ~ GENE)) %>% 
  replace_na(list( GENE = 'NO GENE'))

muts %>% write_csv(here('exploration','mutations_summary_wider.csv'))



#  effects on M and P variants --------------------------------------------


mp_muts = muts %>%
  filter(Strain %in% c('M_strains', 'P_strains'))


mp_muts %>% 
  # remove synonymous variants, not informative
  filter(Variant != 'synonymous_variant') %>%
  group_by(Strain) %>% 
  count(GENE) %>% 
  ggplot(aes(x = fct_reorder(GENE,n), y = n, fill = GENE)) +
  geom_histogram(stat="identity")  +
  labs(x = 'Gene variant',
       y = 'Number of total elements') +
  facet_wrap(~Strain,
             ncol = 1) +
  geom_text(aes(y = n+(2), x = GENE, label = round(n,0))) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(here('exploration','mutations_MP_strains.pdf'), height = 8, width = 10)


mp_muts %>% 
  # remove synonymous variants, not informative
  filter(Variant != 'synonymous_variant') %>%
  group_by(Strain, Colony, GENE) %>%
  summarise(N = n()) %>% 
  ggplot(aes(x = fct_reorder(GENE, N), y = N, fill = GENE)) +
  geom_histogram(stat="identity")  +
  labs(x = 'Gene variant',
       y = 'Number of total elements') +
  facet_wrap(~Strain*Colony)+
  geom_text(aes(y = N+(1), x = GENE, label = round(N,0))) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(here('exploration','mutations_MP_strains_big.pdf'), height = 18, width = 25)

mp_muts %>% 
  distinct(Colony, Strain) %>% 
  group_by(Strain) %>% 
  summary(N = n())

# calculate proportion of mutations per colony
# is a quick way to get a sense of how often they appear
mp_muts %>% 
  # remove synonymous variants, not informative
  filter(Variant != 'synonymous_variant') %>%
  group_by(Strain, GENE) %>%
  summarise(N = n()) %>% 
  group_by(Strain) %>% 
  # 19 P colonies, 12 M colonies
  mutate(prop = case_when(Strain == 'P_strains' ~ N / 19,
                          Strain == 'M_strains' ~ N / 12)) %>% 
  # remove rubish
  # filter(!(GENE %in% c('NO GENE','insB6_1', 'insA6_1'))) %>% 
  ggplot(aes(x = fct_reorder(GENE, prop), y = prop, fill = GENE)) +
  geom_histogram(stat="identity")  +
  labs(x = 'Gene variant',
       y = 'Proportion per colony') +
  geom_text(aes(y = prop+(0.1), x = GENE, label = round(prop,2))) + 
  geom_hline(yintercept = 1) +
  facet_wrap(~Strain,
             ncol = 1)+
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(here('exploration','mut_props_MP_strains.pdf'), height = 8, width = 10)




#  effects on pangenome strains --------------------------------------------


pan_muts = muts %>%
  filter(!(Strain %in% c('M_strains', 'P_strains', 'Nissle')))



pan_muts %>% 
  # remove synonymous variants, not informative
  filter(Variant != 'synonymous_variant') %>%
  group_by(Strain) %>% 
  count(GENE) %>% 
  ggplot(aes(x = fct_reorder(GENE,n), y = n, fill = GENE)) +
  geom_histogram(stat="identity")  +
  labs(x = 'Gene variant',
       y = 'Number of total elements') +
  facet_wrap(~Strain,
             ncol = 3,
             scales = 'free_x') +
  geom_text(aes(y = n+(2), x = GENE, label = round(n,0))) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(here('exploration','mutations_Pangenome_strains.pdf'), height = 8, width = 10)



pan_muts %>% 
  # remove synonymous variants, not informative
  filter(Variant != 'synonymous_variant') %>%
  # group_by(Strain) %>% 
  count(GENE) %>% 
  ggplot(aes(x = fct_reorder(GENE,n), y = n, fill = GENE)) +
  geom_histogram(stat="identity")  +
  labs(x = 'Gene variant',
       y = 'Number of total elements') +
  # facet_wrap(~Strain,
  #            ncol = 3,
  #            scales = 'free_x') +
  geom_text(aes(y = n+(1), x = GENE, label = round(n,0))) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(here('exploration','mutations_Pangenome_total.pdf'), height = 8, width = 11)














