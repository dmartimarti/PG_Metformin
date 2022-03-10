library(tidyverse)
library(ggrepel)
library(here)
library(readr)

gene_hits = read_tsv(here("results", "gene_hits.txt"))


gene_hits %>% 
  filter(avg_beta > 0.8 | maxp > 7) %>% 
  ggplot(aes(x=avg_beta, y=maxp, colour=avg_maf, size=hits, label=gene)) +
  geom_point(alpha=0.5) +
  geom_text_repel(aes(size=30), 
                  show.legend = FALSE, 
                  colour='black',
                  max.overlaps = Inf) +
  scale_size("Number of \nk-mers", range=c(1,10)) +
  scale_colour_gradient('Average \nMAF') +
  theme_bw(base_size=14) +
  xlab("Average effect size") +
  ylab("Maximum -log10(p-value)")


ggsave('worm_phenotype_ALL_no_biofilm.pdf',
       height = 10, width = 13)
