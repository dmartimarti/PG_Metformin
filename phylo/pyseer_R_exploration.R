
# libraries ---------------------------------------------------------------


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





# reading pyseer output ---------------------------------------------------

#   eggnog annotation    #

eggnog = read_csv("eggnog_annotation_fixed.csv") %>% 
  rename(variant = Gene)


worm_all_pyseer = read_delim("results/worm_phenotype_ALL_no_biofilm.tsv", 
           delim = "\t", escape_double = FALSE, 
           trim_ws = TRUE)


# transforming the data to get a full vector of genes

worm_all_pyseer %>% 
  left_join(eggnog) %>% 
  filter(`lrt-pvalue` < 0.01) %>%
  # filter(`filter-pvalue` < 0.00001) %>%
  # filter(str_detect(variant, 'group'))
  mutate(variant = case_when(str_detect(variant, 'group') & Preferred_name != '-' ~ Preferred_name,
                             TRUE ~ variant)) %>% 
  separate_rows(variant, sep = '~~~') %>% 
  mutate(variant = str_replace(variant, pattern = '_[:digit:]', '')) %>% 
  filter(!(str_detect(variant, 'group'))) %>% 
  distinct(variant)  %>% 
  write_delim("worm_ALL.txt")


worm_all_pyseer %>% 
  left_join(eggnog) %>% 
  filter(`lrt-pvalue` < 0.01) %>%
  separate_rows(COG_category, sep = ',') %>% 
  count(COG_category) %>% 
  filter(!(COG_category %in% c('Not annotated', 'S'))) %>% 
  ggplot(aes( x = fct_reorder(COG_category, n), y = n )) +
  geom_col(position = 'dodge')

worm_all_pyseer %>% 
  left_join(eggnog) %>% 
  filter(`lrt-pvalue` < 0.01) %>% 
  separate_rows(KEGG_Pathway, sep= ',') %>% 
  count(KEGG_Pathway) %>% 
  filter(KEGG_Pathway != '-') %>% 
  ggplot(aes( y = fct_reorder(KEGG_Pathway, n), x = n )) +
  geom_col(position = 'dodge')







