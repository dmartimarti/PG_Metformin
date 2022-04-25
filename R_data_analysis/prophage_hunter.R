
# libraries ----------------------------------------------------------------

library(tidyverse)
library(here)
library(readr)
library(readxl)
library(cowplot)
library(glue)

theme_set(theme_cowplot(17))



# load the data -----------------------------------------------------------

metadata = read_excel('metadata.xlsx', sheet = 'selected')

# which genomes to iterate
genomes = metadata %>% filter(Prophage == 'Yes') %>% pull(ID)

# loop over the folders with prophages and filter things below score 0.5
proph = tibble()
for (genome in genomes) {
  temp = read_delim(glue("results_web/{genome}/01.Main_output.txt"), 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE, 
                    col_types = cols()) %>% 
    filter(Score > 0.5) %>% 
    mutate(ID = genome, .before = `Candidate ID`)
  
  proph = proph %>% 
    bind_rows(temp)
}

# how many active prophages do we find
active = proph %>% 
  filter(Category == 'Active') %>% 
  group_by(ID) %>% 
  count() %>% 
  full_join(metadata) %>% 
  mutate(n = replace_na(n, 0)) 

# total prophages
total = proph %>% 
  # filter(Category == 'Active') %>% 
  group_by(ID) %>% 
  count() %>% 
  full_join(metadata) %>% 
  mutate(n = replace_na(n, 0)) 


active %>% 
  filter(Origin == 'ECOREF') %>% 
  ggplot(aes(x = Phenotype, y = n, fill = Phenotype)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.height = 0.2)) + 
  labs(
    x = 'Bacterial phenotype',
    y = 'Prophage in genome'
  )

ggsave(here('exploration','active_prophages_phenotype_ECOREF.pdf'), 
       height = 8, width = 10)

total %>% 
  filter(Origin == 'ECOREF') %>% 
  ggplot(aes(x = Phenotype, y = n, fill = Phenotype)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.height = 0.2)) +
  labs(
    x = 'Bacterial phenotype',
    y = 'Prophage in genome'
  )

ggsave(here('exploration','total_prophages_phenotype_ECOREF.pdf'), 
       height = 8, width = 10)


proph %>% 
  group_by(`Closest phage`) %>% 
  count() %>% 
  filter(`Closest phage` != 'N/A') %>% 
  ggplot(aes(y = fct_reorder(`Closest phage`,n), x = n)) +
  geom_bar(stat = 'identity', fill = 'grey50', color = 'black') +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1, 
                               size = 15)
  ) +
  labs(
    y = 'Closest phage',
    x = 'Number of prophages found'
  )

ggsave(here('exploration','most_common_prophage.pdf'), 
       height = 9, width = 11)


metadata %>% 
  group_by(Prophage) %>% 
  count() %>% 
  ungroup %>% 
  mutate(per = round(100*n/sum(n),2)) %>% 
  ggplot(aes(x = Prophage, y = per)) +
  geom_bar(stat = 'identity', color = 'black', aes(fill = Prophage)) +
  scale_fill_manual(values = c('#F05B22','#1FCCBB')) +
  geom_text(aes(label = paste0(per,'%')), y = 40, size = 6) + 
  guides(fill = 'none')

ggsave(here('exploration','barplot_prophage_presence.pdf'), 
       height = 8, width = 7)



proph %>% 
  # filter(Score > 0.8) %>% 
  full_join(metadata) %>% 
  group_by(`Closest phage`, Phenotype) %>% 
  summarise(N = n()) %>% 
  ungroup %>% 
  group_by(`Closest phage`) %>% 
  mutate(prop = N/sum(N)) %>% 
  filter(`Closest phage` != 'N/A') %>% 
  drop_na(`Closest phage`) %>% 
  ggplot(aes(x = prop, y = `Closest phage`, fill = Phenotype)) +
  geom_bar(stat = 'identity')
  # facet_wrap(~Phenotype)

ggsave(here('exploration','barplot_prophage_presence_phenotype.pdf'), 
       height = 8, width = 13)
