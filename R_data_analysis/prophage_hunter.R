
# libraries ----------------------------------------------------------------

library(tidyverse)
library(here)
library(readr)
library(readxl)
library(cowplot)

theme_set(theme_cowplot(17))



# load the data -----------------------------------------------------------

metadata = read_excel('metadata.xlsx', sheet = 'selected')

# which genomes to iterate
genomes = metadata %>% filter(Prophage == 'Yes') %>% pull(ID)

# loop over the folders with prophages and filter things below score 0.5
proph = tibble()
for (genome in genomes) {
  temp = read_delim("results_web/NT12001/01.Main_output.txt", 
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
  ggplot(aes(x = `Closest phage`, y = n)) +
  geom_bar(stat = 'identity') +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, 
                               size = 8)
  ) +
  labs(
    x = 'Closes phage',
    y = 'Number of prophages found'
  )

ggsave(here('exploration','most_common_prophage.pdf'), 
       height = 8, width = 10)


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






