
# libraries ---------------------------------------------------------------

library(tidyverse)
library(here)
library(readr)
library(readxl)
library(cowplot)
library(ggpubr)
library(showtext)
library(ggtext)


theme_set(theme_cowplot(15))



theme_nice_45 = theme_set(theme_cowplot(16) + 
                            theme(
                              axis.text.x = element_text(angle = 45, vjust = 0.5),
                              plot.title = element_textbox_simple(family = 'patua-one', size = 20),
                              plot.title.position = 'plot',
                              plot.caption = element_markdown(hjust = 0, color='grey50',
                                                              margin = margin(t=10)),
                              plot.caption.position = 'plot'
                            ))


# load data ---------------------------------------------------------------


metadata = read_excel("~/Documents/MRC_postdoc/Pangenomic/metadata/MAIN_metadata.xlsx", 
                            sheet = "metadata")

prism = read_csv("summary.csv") %>% 
  mutate(genome = as.factor(genome)) 

prism_smiles = read_csv("summary_smiles.csv") %>% 
  mutate(genome = as.factor(genome)) 


genome_info = read_delim("~/Documents/MRC_Postdoc/Pangenomic/pangenome_analysis/ALL/phylo_analysis/assemblies/no_evo/quast_quality/transposed_report.tsv", 
                         delim = "\t", escape_double = FALSE, 
                         col_types = cols(Assembly = col_character()), 
                         trim_ws = TRUE) %>% 
  rename(genome = Assembly,
         genome_length = `Total length (>= 0 bp)`,
         genome_length_500 =`Total length`,
         gc_content = `GC (%)`) %>% 
  select(genome, genome_length, genome_length_500, gc_content)


# exploration plots -------------------------------------------------------


prism %>% 
  group_by(genome) %>% 
  count() %>% 
  ggplot(aes(n)) +
  geom_histogram(aes(y = ..density..) ,bins = 8, color = 'black',
                 fill = 'grey50') +
  # geom_density( fill = 'darkslategray3', alpha = .5) +
  labs(x = 'Number of BGCs',
       y = 'Density') +
  theme_cowplot(19)

ggsave(here('exploration', 'BGC_density.pdf'), height = 8, width = 10)


### number of BGCs ####

prism %>% 
  count(type) %>% 
  ggplot(aes(x = fct_reorder(type, n,.desc = TRUE), y = n, fill = type)) +
  geom_bar(stat='identity', color = 'black') +
  theme_cowplot(14) +
  labs(y = 'Number of BGCs',
       x = NULL, 
       # title = 'Distribution of BGCs',
       ) +
  guides(fill = 'none') +
  theme_cowplot(16) +
  theme_nice_45

ggsave(here('exploration', 'BGC_numbers.pdf'), height = 8, width = 10)



prism %>% 
  separate_rows(type, sep = '\\|') %>% 
  count(type) %>% 
  ggplot(aes(x = fct_reorder(type, n,.desc = TRUE), y = n, fill = type)) +
  geom_bar(stat='identity', color = 'black') +
  theme_cowplot(16) +
  labs(y = 'Number of BGCs',
       x = NULL) +
  guides(fill = 'none') +
  # viridis::scale_fill_viridis(discrete = T) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))


ggsave(here('exploration', 'BGC_separate_numbers.pdf'), height = 8, width = 10)


## BGCs per genome size ####

prism %>% 
  group_by(genome) %>% 
  count() %>% 
  left_join(genome_info) %>% 
  filter(!(genome %in% c('2.3', '6'))) %>% 
  mutate(genome_length = genome_length/1000000) %>% 
  ggplot(aes(x = genome_length, y = n)) +
  geom_point(shape = 21, fill = 'dodgerblue', size = 3, alpha = 0.5) +
  geom_smooth(method = 'lm') +
  labs(x = 'Genome size (Mb)',
       y = 'Number of BGCs') +
  # scale_x_continuous(limits = c(4.3,6.516906)) +
  # scale_y_continuous(limits = c(1,9), breaks = c(2,4,6,8)) +
  stat_cor(method = "pearson", label.x = 5, label.y = 9,
           p.accuracy = 0.001, r.accuracy = 0.01, size = 7,
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_cowplot(17)

ggsave(here('exploration', 'genome_size_vs_BGCs.pdf'), height = 8, width = 10)

## GC vs number of BGCs ####

library(ggpubr)

prism %>% 
  group_by(genome) %>% 
  count() %>% 
  left_join(genome_info) %>% 
  filter(!(genome %in% c('2.3', '6')))  %>% 
  ggplot(aes(x = gc_content, y = n)) +
  geom_point(shape = 21, fill = 'dodgerblue', size = 3, alpha = 0.5) +
  geom_smooth(method = 'lm') +
  labs(x = 'GC content (%)',
       y = 'Number of BGCs') +
  scale_y_continuous(breaks = c(2,4,6,8)) +
  stat_cor(method = "pearson", label.x = 50.4, label.y = 9,
           p.accuracy = 0.001, r.accuracy = 0.01, size = 7,
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_cowplot(17)

ggsave(here('exploration', 'gc_content_vs_BGCs.pdf'), height = 8, width = 10)



### repeated elements within genomes ####

prism %>% 
  # separate_rows(type, sep = '\\|') %>%
  group_by(genome, type) %>% 
  count()  %>% 
  mutate(repeated = ifelse(n > 1, 'yes', 'no')) %>% 
  group_by(type, repeated) %>% 
  count() %>% 
  filter(repeated == 'yes') %>% 
  ggplot(aes(x = fct_reorder(type, n, .desc=TRUE), y = n, 
             fill = repeated, group = repeated)) +
  geom_bar(stat='identity', fill = 'dodgerblue3', color = 'black') +
  labs(
    x = 'BGC type',
    y = 'Number of times repeated within a genome'
  ) +
  theme_cowplot(13)  +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

ggsave(here('exploration', 'BGC_repeated.pdf'), height = 5, width = 6)


## smiles ####

smiles.sum = prism_smiles %>% 
  select(genome:mol_weight) %>% 
  drop_na(smiles) %>% 
  group_by(genome, cluster) %>% 
  distinct(smiles, .keep_all = T) %>% 
  ungroup %>% 
  count(smiles) %>% 
  arrange(desc(n))

smiles.sum = smiles.sum %>% 
  left_join(prism_smiles %>% select(smiles, mol_weight)) %>% 
  distinct(smiles, .keep_all = T)

smiles.sum %>% 
  write_csv(here('exploration', 'smiles_summary.csv'))

smiles.sum %>% 
  arrange(desc(n)) %>% 
  mutate(num_name = seq(1,dim(smiles.sum)[1],1)) %>%
  filter(num_name < 30) %>% 
  mutate(num_name = as.factor(num_name)) %>% 
  ggplot(aes(x = fct_reorder(num_name, n, .desc = T), y = n)) +
  geom_bar(stat='identity') +
  theme_cowplot(17) +
  labs(y = 'Compound count',
       x = 'Compound  (simplified)') 

ggsave(here('exploration', 'smiles_numbers.pdf'), height = 8, width = 10)


smiles.sum %>% 
  ggplot(aes(mol_weight)) +
  geom_density()

smiles.sum %>% 
  arrange(desc(n)) %>% 
  mutate(num_name = seq(1,dim(smiles.sum)[1],1)) %>% 
  write_csv(here('exploration', 'smiles_summary.csv'))

