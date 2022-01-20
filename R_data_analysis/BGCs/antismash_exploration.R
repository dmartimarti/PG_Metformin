
# libraries ---------------------------------------------------------------


library(tidyverse)
library(readxl)
library(readr)
library(cowplot)
library(here)
library(ggpubr)
library(showtext)
library(ggtext)

# Theme settings for all plots

theme_nice = theme_set(theme_cowplot(16) + 
            theme(
              plot.title = element_textbox_simple(family = 'patua-one', size = 20),
              plot.title.position = 'plot',
              plot.caption = element_markdown(hjust = 0, color='grey50',
                                              margin = margin(t=10)),
              plot.caption.position = 'plot'
            ))

theme_nice_45 = theme_set(theme_cowplot(16) + 
                         theme(
                           axis.text.x = element_text(angle = 45, vjust = 0.5),
                           plot.title = element_textbox_simple(family = 'patua-one', size = 20),
                           plot.title.position = 'plot',
                           plot.caption = element_markdown(hjust = 0, color='grey50',
                                                           margin = margin(t=10)),
                           plot.caption.position = 'plot'
                         ))



font_add_google('Patua One', 'patua-one')

showtext_auto()

# load data ---------------------------------------------------------------


bgc = read_csv("antismash_summary.csv")


smiles = read_csv("antismash_summary_smiles.csv")


metadata = read_excel("D:/MRC_Postdoc/Pangenomic/pangenome_analysis/metadata/MAIN_metadata.xlsx", 
                      sheet = "metadata")


genome_info = read_delim("D:/MRC_Postdoc/Pangenomic/pangenome_analysis/ALL/phylo_analysis/assemblies/no_evo/quast_quality/transposed_report.tsv", 
                         delim = "\t", escape_double = FALSE, 
                         col_types = cols(Assembly = col_character()), 
                         trim_ws = TRUE) %>% 
  rename(genome = Assembly,
         genome_length = `Total length (>= 0 bp)`,
         genome_length_500 =`Total length`,
         gc_content = `GC (%)`) %>% 
  select(genome, genome_length, genome_length_500, gc_content)



# this is the info of the core genome genes
core_genome_info = read_delim("D:/MRC_Postdoc/Pangenomic/pangenome_analysis/ALL/phylo_analysis/panaroo_results_noEVO/core_genomes/quast_quality/transposed_report_core_genomes.tsv", 
                         delim = "\t", escape_double = FALSE, 
                         col_types = cols(Assembly = col_character()), 
                         trim_ws = TRUE) %>% 
  rename(genome = Assembly,
         genome_length = `Total length (>= 0 bp)`,
         genome_length_500 =`Total length`,
         gc_content_core = `GC (%)`) %>% 
  select(genome, genome_length, genome_length_500, gc_content_core)


# clean up the summary file

bgc %>% 
  filter(genome == 98)

# fix the order of some names
bgc = bgc %>% 
  mutate(type = case_when(type == 'T1PKS|NRPS' ~ 'NRPS|T1PKS',
                          type == 'NRPS-like|NRPS' ~ 'NRPS|NRPS-like',
                          type == 'ladderane|arylpolyene' ~ 'arylpolyene|ladderane',
                          type == 'NRPS|T1PKS|NRPS-like' ~ 'NRPS|NRPS-like|T1PKS',
                          type == 'T1PKS|NRPS|NRPS-like' ~ 'NRPS|NRPS-like|T1PKS',
                          TRUE ~ type))



# exploration  ############

### number of BGCs ####

bgc %>% 
  count(type) %>% 
  ggplot(aes(x = fct_reorder(type, n,.desc = TRUE), y = n, fill = type)) +
  geom_bar(stat='identity', color = 'black') +
  theme_cowplot(14) +
  labs(y = 'Number of BGCs',
       x = NULL, 
       # title = 'Distribution of BGCs',
       caption = '<i> Distribution of the different BGCs across the 746 
       bacterial strains in our collection </i>') +
  guides(fill = 'none') +
  theme_cowplot(16) +
  theme_nice 

ggsave(here('exploration', 'BGC_numbers.pdf'), height = 8, width = 10)



bgc %>% 
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


## BGC density ####
bgc %>% 
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



##  genome size ####
genome_info %>% 
  filter(genome != '2.3') %>% 
  mutate(genome_length = genome_length/1000000) %>%
  ggplot(aes(genome_length)) +
  geom_histogram(aes(y = ..density..) ,bins = 30, color = 'black',
                 fill = 'firebrick2') +
  geom_density(fill = 'darkslategray3', alpha = .3) +
  labs(x = 'Genome size (Mb)',
       y = 'Density') +
  theme_cowplot(14)

ggsave(here('exploration', 'genome_size_density.pdf'), height = 6, width = 8)

# info about the core genome
core_genome_info %>% 
  filter(genome != '2.3') %>% 
  mutate(genome_length = genome_length/1000000) %>%
  ggplot(aes(genome_length)) +
  geom_histogram(aes(y = ..density..) ,bins = 30, color = 'black',
                 fill = 'firebrick2') +
  geom_density(fill = 'darkslategray3', alpha = .3) +
  labs(x = 'Genome size (Mb)',
       y = 'Density') +
  theme_cowplot(14)

ggsave(here('exploration', 'genome_size_CORE_density.pdf'), height = 6, width = 8)

##  GC content ####
genome_info %>% 
  filter(genome != '2.3') %>% 
  ggplot(aes(gc_content)) +
  geom_histogram(aes(y = ..density..) ,bins = 30, color = 'black',
                 fill = c("#228B22")) +
  geom_density(fill = c("#228B22"), alpha = .3) +
  labs(x = 'GC content (%)',
       y = 'Density') +
  theme_cowplot(17)

ggsave(here('exploration', 'GC_content_density.pdf'), height = 8, width = 10)

core_genome_info %>% 
  filter(genome != '2.3') %>% 
  ggplot(aes(gc_content_core)) +
  geom_histogram(aes(y = ..density..) ,bins = 30, color = 'black',
                 fill = c("#228B22")) +
  geom_density(fill = c("#228B22"), alpha = .3) +
  labs(x = 'GC content (%)',
       y = 'Density') +
  theme_cowplot(17)


ggsave(here('exploration', 'GC_content_CORE_density.pdf'), height = 8, width = 10)

##  GC vs genome size ####
genome_info %>% 
  filter(genome != '2.3') %>% 
  mutate(genome_length = genome_length/1000000) %>%
  ggplot(aes(x = genome_length, y = gc_content)) +
  geom_point(shape = 21, fill = 'dodgerblue', size = 3, alpha = 0.5) +
  geom_smooth(method = 'lm') +
  labs(x = 'Genome size (Mb)',
       y = 'GC content (%)') +
  # scale_x_continuous(limits = c(4.3,6.516906)) +
  # scale_y_continuous(limits = c(1,9), breaks = c(2,4,6,8)) +
  stat_cor(method = "pearson", label.x = 5.2, label.y = 51.1,
           p.accuracy = 0.001, r.accuracy = 0.01, size = 7,
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_cowplot(17)

ggsave(here('exploration', 'GC_vs_genomeSize.pdf'), height = 8, width = 10)


genome_info %>% 
  left_join(core_genome_info %>% 
              select(genome, gc_content_core)) %>% 
  filter(genome != '2.3') %>% 
  mutate(genome_length = genome_length/1000000) %>%
  ggplot(aes(x = genome_length, y = gc_content_core)) +
  geom_point(shape = 21, fill = 'dodgerblue', size = 3, alpha = 0.5) +
  geom_smooth(method = 'lm') +
  labs(x = 'Genome size (Mb)',
       y = 'GC content (%)') +
  # scale_x_continuous(limits = c(4.3,6.516906)) +
  # scale_y_continuous(limits = c(1,9), breaks = c(2,4,6,8)) +
  stat_cor(method = "pearson", label.x = 5.2, label.y = 52.7,
           p.accuracy = 0.001, r.accuracy = 0.01, size = 7,
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_cowplot(17)

ggsave(here('exploration', 'GC_vs_genomeSize_CORE.pdf'), height = 8, width = 10)


## genome size vs number of BGCs ####

bgc %>% 
  group_by(genome) %>% 
  count() %>% 
  left_join(genome_info) %>% 
  filter(genome != '2.3') %>% 
  mutate(genome_length = genome_length/1000000) %>% 
  ggplot(aes(x = genome_length, y = n)) +
  geom_point(shape = 21, fill = 'dodgerblue', size = 3, alpha = 0.5) +
  geom_smooth(method = 'lm') +
  labs(x = 'Genome size (Mb)',
       y = 'Number of BGCs') +
  scale_x_continuous(limits = c(4.3,6.516906)) +
  scale_y_continuous(limits = c(1,9), breaks = c(2,4,6,8)) +
  stat_cor(method = "pearson", label.x = 5, label.y = 9,
           p.accuracy = 0.001, r.accuracy = 0.01, size = 7,
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) +
  theme_cowplot(17)

ggsave(here('exploration', 'genome_size_vs_BGCs.pdf'), height = 8, width = 10)


## GC vs number of BGCs ####

library(ggpubr)

bgc %>% 
  group_by(genome) %>% 
  count() %>% 
  left_join(genome_info) %>% 
  filter(genome != '2.3') %>% 
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

bgc %>% 
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

smiles.sum = smiles %>% 
  select(genome:mol_weight) %>% 
  drop_na(smiles) %>% 
  group_by(genome, cluster) %>% 
  distinct(smiles, .keep_all = T) %>% 
  ungroup %>% 
  count(smiles) %>% 
  arrange(desc(n))

smiles.sum = smiles.sum %>% 
  left_join(smiles %>% select(smiles, mol_weight)) %>% 
  distinct(smiles, .keep_all = T)

smiles.sum %>% 
  write_csv(here('exploration', 'smiles_summary.csv'))

smiles.sum %>% 
  arrange(desc(n)) %>% 
  mutate(num_name = seq(1,28,1),
         num_name = as.factor(num_name)) %>% 
  ggplot(aes(x = fct_reorder(num_name, n, .desc = T), y = n)) +
  geom_bar(stat='identity') +
  theme_cowplot(17) +
  labs(y = 'Compound count',
       x = 'Compound  (simplified)') 

ggsave(here('exploration', 'smiles_numbers.pdf'), height = 8, width = 10)



# add metadata ------------------------------------------------------------


bgc.meta = bgc %>% 
  left_join(metadata %>% 
              mutate(genome = str_sub(fasta,1, -7)) %>% 
              select(genome, phylogroup, Broadphenotype))

bgc.meta %>% 
  separate_rows(type, sep = '\\|') %>% 
  group_by(phylogroup) %>% 
  count(type) %>% 
  ggplot(aes(x = fct_reorder(type, n, .desc=TRUE), y = n, fill = type)) +
  geom_bar(stat='identity', color = 'black') +
  theme_cowplot(12) +
  panel_border() +
  labs(y = 'Number of elements',
       x = NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  facet_wrap(~phylogroup, ncol = 2)

ggsave(here('exploration', 'BGC_by_phylogroup.pdf'), height = 8, width = 6)


bgc.meta %>% 
  filter(phylogroup != 'cladeI') %>% 
  separate_rows(type, sep = '\\|') %>% 
  group_by(genome, phylogroup, type) %>% 
  summarise(N = n()) %>% 
  group_by(phylogroup) %>% 
  mutate(
    new_n = n(),
    prop = N / n()) %>% 
  group_by(phylogroup, type) %>% 
  summarise(N = sum(prop)) %>% 
  ggplot(aes(x = fct_reorder(type, N, .desc=T), y = N, fill = type)) +
  geom_bar(stat='identity', color = 'black') +
  theme_cowplot(16) +
  panel_border() +
  labs(y = 'Proportion of elements',
       x = NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"))) +
  facet_wrap(~phylogroup, ncol = 1) +
  guides(fill = 'none')


ggsave(here('exploration', 'BGC_prop_by_phylogroup.pdf'), height = 12, width = 6)







# co-ocurrence ------------------------------------------------------------

library(cooccur)
library(visNetwork)


### create my own presence/absence matrix to analyse
# genomes as rows, bgcs as columns

bgc_pa = bgc %>% 
  separate_rows(type, sep = '\\|') %>% 
  select(genome, type) %>% 
  mutate(present = 1) %>% 
  pivot_wider(names_from = type, values_from = present, 
              values_fn = {mean}, values_fill = 0)

genome_names = bgc_pa %>% pull(genome)

bgc_pa_matrix = bgc_pa %>% select(-genome) %>% 
  as.matrix
   
rownames(bgc_pa_matrix) = genome_names



### co-ocurrence analysis

bgc_pa_matrix = t(bgc_pa_matrix)

co_bgc = print(cooccur(bgc_pa_matrix, spp_names = TRUE))

# Check sp1_name matches numeric label for species.
co_bgc[, 'sp1_name'] == rownames(bgc_pa_matrix)[co_bgc$sp1]
co_bgc[, 'sp2_name'] == rownames(bgc_pa_matrix)[co_bgc$sp2]


# Create a data frame of the nodes in the network. 
nodes = data.frame(id = 1:nrow(bgc_pa_matrix),
                    label = rownames(bgc_pa_matrix),
                    color = "#606482",
                    shadow = TRUE)

# Create an edges dataframe from the significant pairwise co-occurrences.
edges = data.frame(from = co_bgc$sp1, to = co_bgc$sp2,
                    color = ifelse(co_bgc$p_lt <= 0.05, "#F2D133", "#0402A6"),
                    dashes = ifelse(co_bgc$p_lt <= 0.05, FALSE, TRUE))

# Plot.
visNetwork(nodes = nodes, edges = edges) %>%
  visIgraphLayout(layout = "layout_with_kk")


write.csv(as.data.frame(bgc_pa_matrix), 'bgc_pa_matrix.csv')
