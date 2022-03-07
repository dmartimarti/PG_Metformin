
# libraries ---------------------------------------------------------------



library(ape)
library(tidyverse)
library(readr)
library(cowplot)
library(here)

theme_set(theme_cowplot(15))

# check if names are correct for input files ------------------------------



tree = ape::read.tree('tree.newick')
my_file = read_delim("biofilm/dbgwas_biofilm_PG.txt", delim = '\t')

tree_labels =  tree$tip.label
my_labels = my_file$ID


setdiff(my_labels, tree_labels)





# results exploration -----------------------------------------------------

## Biofilm ####

annotations = read_delim("results/biofilms/textualOutput/all_comps_annotations_info.tsv", 
                                         delim = "\t", escape_double = FALSE, 
                                         trim_ws = TRUE)

nodes = read_delim("results/biofilms/textualOutput/all_comps_nodes_info.tsv", 
                                        delim = "\t", escape_double = FALSE, 
                                        trim_ws = TRUE)



nodes %>% 
  filter(`Significant?` == 'Yes')



nodes %>% 
  ggplot(aes(CompId)) + 
  geom_bar(aes(fill = `Significant?`)) + 
  labs(
    y = 'Number of nodes',
    x = 'Component ID',
    caption = 'Number of nodes per component from the DBGWAS output'
    )

ggsave(here('results', 'biofilms', 'R_exploration', 'nodes_count.pdf'),
       height = 9, width = 12)


nodes %>% 
  group_by(CompId) %>% 
  count() %>% 
  arrange(desc(n)) %>%
  ggplot(aes(x = n)) +
  geom_histogram(color = 'black', fill = 'grey70') +
  labs(
    y = 'Node count',
    x = 'Components',
    caption = 'Node count histogram'
  )

ggsave(here('results', 'biofilms', 'R_exploration', 'nodes_histogram.pdf'),
       height = 8, width = 11)



nodes %>% 
  filter(`Significant?` == 'Yes') %>% 
  arrange(EstEffect) %>%
  mutate(Effect = case_when(EstEffect > 0 ~ 'Biofilm',
                               EstEffect < 0 ~ 'WT')) %>% 
  ggplot(aes(x = EstEffect, 
             y = fct_reorder(NodeId, EstEffect),
             fill = Effect)) +
  geom_bar(stat = 'identity') +
  labs(
    y = 'Node ID',
    x = 'Estimated effect'
  )

ggsave(here('results', 'biofilms', 'R_exploration', 'EstEffect_sig.pdf'),
       height = 22, width = 9)



nodes %>% 
  filter(`Significant?` == 'Yes') %>% 
  arrange(EstEffect) %>%
  drop_na(`Annotations(sep=~~~)`) %>% 
  mutate(Effect = case_when(EstEffect > 0 ~ 'Biofilm',
                            EstEffect < 0 ~ 'WT')) %>% 
  ggplot(aes(x = EstEffect, 
             y = fct_reorder(NodeId, EstEffect),
             fill = Effect)) +
  geom_bar(stat = 'identity') +
  labs(
    y = 'Node ID',
    x = 'Estimated effect'
  )

ggsave(here('results', 'biofilms', 'R_exploration', 'EstEffect_sig_annotated.pdf'),
       height = 15, width = 9)

