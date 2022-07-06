library(tidyverse)
library(ggrepel)
library(here)
library(readr)
library(cowplot)

gene_hits = read_tsv(here("results", "gene_hits.txt"))


gene_hits %>% 
  mutate(gene_name = case_when(avg_beta > 1 | 
                                 maxp > 8 ~ gene,
                               str_detect(gene, 'mdt') ~ gene)) %>% 
  # mutate(label_size = case_when(gene_name != 'mdtB' ~ 20,
  #                               gene_name == 'mdtB' ~ 65)) %>% 
  mutate(label_size = case_when(str_detect(gene_name, 'mdt')  ~ 65,
                                TRUE ~ 20)) %>% 
  # filter(avg_beta > 0.8 | maxp > 7) %>% 
  ggplot(aes(x=avg_beta, y=maxp, colour=avg_maf, size=hits, label=gene_name)) +
  geom_point(alpha=0.5) +
  geom_text_repel(aes(size=label_size), 
                  show.legend = FALSE, 
                  colour='black',
                  max.overlaps = Inf) +
  scale_size("Number of \nk-mers", range=c(1,10)) +
  scale_colour_gradient('Average \nMAF') +
  theme_cowplot(15) +
  xlab("Average effect size") +
  ylab("Maximum -log10(p-value)")


ggsave('worm_phenotype_ALL_no_biofilm.pdf',
       height = 6, width = 8)




gene_hits %>% 
  mutate(gene_name = case_when(avg_beta > 1 | 
                                 maxp > 8 ~ gene,
                               str_detect(gene, 'mdt') ~ gene)) %>% 
  # mutate(label_size = case_when(gene_name != 'mdtB' ~ 20,
  #                               gene_name == 'mdtB' ~ 65)) %>% 
  mutate(label_size = case_when(str_detect(gene_name, 'mdt')  ~ 65,
                                TRUE ~ 20)) %>% 
  # filter(avg_beta > 0.8 | maxp > 7) %>% 
  ggplot(aes(x=avg_beta, y=maxp, colour=avg_maf, size=hits, label=gene_name)) +
  geom_point(alpha=0.5) +
  geom_text_repel(aes(size=label_size), 
                  show.legend = FALSE, 
                  colour='black',
                  max.overlaps = Inf) +
  scale_size("Number of \nk-mers", range=c(1,10)) +
  scale_colour_gradient('Average \nMAF') +
  theme_cowplot(15) +
  xlab("Average effect size") +
  ylab("Maximum -log10(p-value)")

ggsave('worm_phenotype_ALL_no_biofilm_v2.pdf',
       height = 6, width = 8)


