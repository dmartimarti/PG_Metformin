
# libraries ---------------------------------------------------------------


library(tidyverse)
library(readxl)
library(readr)
library(cowplot)
library(here)
library(ggpubr)


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

# exploration  ############

bgc %>% 
  count(type) %>% 
  ggplot(aes(x = fct_reorder(type, n,.desc = TRUE), y = n)) +
  geom_bar(stat='identity') +
  theme_cowplot(12) +
  labs(y = 'Number of BGCs',
       x = NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ggsave(here('exploration', 'BGC_numbers.pdf'), height = 8, width = 10)


## BGC density ####
bgc %>% 
  group_by(genome) %>% 
  count() %>% 
  ggplot(aes(n)) +
  geom_histogram(aes(y = ..density..) ,bins = 8, color = 'black',
                 fill = 'dodgerblue2') +
  # geom_density( fill = 'darkslategray3', alpha = .5) +
  labs(x = 'Number of BGCs',
       y = 'Density') +
  theme_cowplot(14)

ggsave(here('exploration', 'BGC_density.pdf'), height = 8, width = 10)



##  genome size ####
genome_info %>% 
  filter(genome != '2.3') %>% 
  ggplot(aes(genome_length)) +
  geom_histogram(aes(y = ..density..) ,bins = 30, color = 'black',
                 fill = 'dodgerblue2') +
  geom_density(fill = 'darkslategray3', alpha = .3) +
  labs(x = 'Genome size',
       y = 'Density') +
  theme_cowplot(14)

ggsave(here('exploration', 'genome_size_density.pdf'), height = 8, width = 10)

# info about the core genome
core_genome_info %>% 
  filter(genome != '2.3') %>% 
  ggplot(aes(genome_length)) +
  geom_histogram(aes(y = ..density..) ,bins = 30, color = 'black',
                 fill = 'dodgerblue2') +
  geom_density(fill = 'darkslategray3', alpha = .3) +
  labs(x = 'Genome size',
       y = 'Density') +
  theme_cowplot(14)


##  GC content ####
genome_info %>% 
  filter(genome != '2.3') %>% 
  ggplot(aes(gc_content)) +
  geom_histogram(aes(y = ..density..) ,bins = 30, color = 'black',
                 fill = '#EB4D4D') +
  geom_density(fill = c("#D64B4B"), alpha = .3) +
  labs(x = 'GC content (%)',
       y = 'Density') +
  theme_cowplot(17)

ggsave(here('exploration', 'GC_content_density.pdf'), height = 8, width = 10)

core_genome_info %>% 
  filter(genome != '2.3') %>% 
  ggplot(aes(gc_content_core)) +
  geom_histogram(aes(y = ..density..) ,bins = 30, color = 'black',
                 fill = '#EB4D4D') +
  geom_density(fill = c("#D64B4B"), alpha = .3) +
  labs(x = 'GC content (%)',
       y = 'Density') +
  theme_cowplot(17)


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



## genome size vs number of BGCs ####

library(ggpubr)

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





## smiles ####

smiles.sum = smiles %>% 
  select(genome:mol_weight) %>% 
  drop_na(smiles) %>% 
  group_by(genome, cluster) %>% 
  distinct(smiles, .keep_all = T) %>% 
  ungroup %>% 
  count(smiles)

smiles.sum %>% 
  ggplot(aes(x = fct_reorder(smiles,n), y = n)) +
  geom_bar(stat='identity') +
  theme_cowplot(12) +
  labs(y = 'Number of times a compound appears',
       x = NULL) +
  theme(
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 

ggsave(here('exploration', 'smiles_numbers.pdf'), height = 8, width = 10)



# add metadata ------------------------------------------------------------


bgc.meta = bgc %>% 
  left_join(metadata %>% 
              mutate(genome = str_sub(fasta,1, -7)) %>% 
              select(genome, phylogroup, Broadphenotype))

bgc.meta %>% 
  group_by(phylogroup) %>% 
  count(type) %>% 
  ggplot(aes(x = fct_reorder(type, n), y = n)) +
  geom_bar(stat='identity') +
  theme_cowplot(12) +
  labs(y = 'Number of elements',
       x = NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  facet_wrap(~phylogroup, ncol = 2)

ggsave(here('exploration', 'BGC_by_phylogroup.pdf'), height = 9, width = 14)


bgc.meta %>% 
  filter(phylogroup != 'cladeI') %>% 
  group_by(genome, phylogroup, type) %>% 
  summarise(N = n()) %>% 
  group_by(phylogroup) %>% 
  mutate(
    new_n = n(),
    prop = N / n()) %>% 
  ggplot(aes(x = fct_reorder(type, prop), y = prop)) +
  geom_bar(stat='identity') +
  theme_cowplot(12) +
  labs(y = 'Proportion of total elements',
       x = NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"))) +
  facet_wrap(~phylogroup, ncol = 1)


ggsave(here('exploration', 'BGC_prop_by_phylogroup.pdf'), height = 13, width = 7)







# co-ocurrence ------------------------------------------------------------

library(cooccur)
library(visNetwork)
# Load finches data set.
data(finches)
finches[1:5, 1:5]


co <- print(cooccur(finches, spp_names = TRUE))

# Check sp1_name matches numeric label for species.
co[, 'sp1_name'] == rownames(finches)[co$sp1]
co[, 'sp2_name'] == rownames(finches)[co$sp2]


# Create a data frame of the nodes in the network. 
nodes <- data.frame(id = 1:nrow(finches),
                    label = rownames(finches),
                    color = "#606482",
                    shadow = TRUE)

# Create an edges dataframe from the significant pairwise co-occurrences.
edges <- data.frame(from = co$sp1, to = co$sp2,
                    color = ifelse(co$p_lt <= 0.05, "#B0B2C1", "#3C3F51"),
                                   dashes = ifelse(co$p_lt <= 0.05, TRUE, FALSE))

# Plot.
visNetwork(nodes = nodes, edges = edges) %>%
  visIgraphLayout(layout = "layout_with_kk")





# gsub(pattern = "\\.|\\:|", replacement = ";", x = "R0.2021-03-17:Col0<F>-")
# 
# 
# pattern = 'R[0-9]\\.\\d{4}\\-\\d{2}\\-\\d{2}\\:[:alnum:]+\\<[:alnum:]+\\>.'
# 
# str_subset(string = "R0.2021-03-17:Col0<F>-", pattern = '^R', negate = F)
# 
# 
# 
# str_replace_all(string = "R0.2021-03-17:Col0<F>-", 
#             pattern = "[^\\w0-9-+]", 
#             replacement = ';') %>% 
#   str_split(pattern = ';', simplify = T)
#   




