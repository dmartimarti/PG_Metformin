

# this scripts prepares the data needed to run the ANN python script

library(tidyverse)
library(readxl)
library(openxlsx)

# metadata
metadata = read_excel("~/Documents/MRC_postdoc/Pangenomic/metadata/MAIN_metadata.xlsx")

# worm phenotype
worm = read_excel("~/Documents/MRC_postdoc/Pangenomic/pangenome_analysis/ALL/worm_imaging/ALL_worm_FC.xlsx")

# gene P/A
gene_pa = read_delim("~/Documents/MRC_postdoc/Pangenomic/pangenome_analysis/ALL/phylo_analysis/panaroo_results_noEVO/gene_presence_absence.Rtab", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)



###
# now filter the relevant genes:
different_genes = pyseer_ALL_no_biofilm$variant
# and get only the genes that are significant in either the two categories
pval_thres = 0.0001
lrt_thres = 0.001
sig_genes = pyseer_ALL_no_biofilm %>% 
  filter((`filter-pvalue` < pval_thres) | (`lrt-pvalue` < lrt_thres)) %>% 
  pull(variant)


# we need to transpose the gene_pa to have genomes in a column and genes as rows

colnames(gene_pa)

gene_pa_all = gene_pa %>% 
  filter(Gene %in% different_genes) %>% 
  pivot_longer(cols = `100`:SPC_4, names_to = 'genome',
               values_to = 'presence') %>% 
  pivot_wider(names_from = 'Gene',
              values_from = presence)


gene_pa_sig = gene_pa %>% 
  filter(Gene %in% sig_genes) %>% 
  pivot_longer(cols = `100`:SPC_4, names_to = 'genome',
               values_to = 'presence') %>% 
  pivot_wider(names_from = 'Gene',
              values_from = presence)

dim(gene_pa_t)




## Merging the worm phenotype, metadata and genetic information together ####



# case for all genes included in pyseer
worm_FC_PA_all = metadata %>% 
  mutate(genome = str_sub(fasta, 1, -7),
         .before= 'ID') %>% 
  filter(genome %in% input_pyseer$IDs) %>% 
  distinct(genome, .keep_all = TRUE) %>% 
  filter(Discard == 'No') %>% 
  left_join(input_pyseer %>% 
              rename(genome = IDs)) %>% 
  left_join(gene_pa_all) %>% 
  mutate(FC_class = ntile(FC_worm, 4),
         .before = 'ID')

# case for only genes that are significant
worm_FC_PA_sig = metadata %>% 
  mutate(genome = str_sub(fasta, 1, -7),
         .before= 'ID') %>% 
  filter(genome %in% input_pyseer$IDs) %>% 
  distinct(genome, .keep_all = TRUE) %>% 
  filter(Discard == 'No') %>% 
  left_join(input_pyseer %>% 
              rename(genome = IDs)) %>% 
  left_join(gene_pa_sig) %>% 
  mutate(FC_class = ntile(FC_worm, 4),
         .before = 'ID')



# divide the data in quantiles
quants = worm_FC_PA_sig %>% 
  summarise(FC_class = scales::percent(c(.20,.80)),
            class = quantile(FC_worm, c(.20,.80)))



# save the full model 
worm_FC_PA_all %>% 
  mutate(FC_class = case_when(FC_worm <= quants$class[1] ~ 'low',
                              FC_worm > quants$class[1] & 
                                FC_worm <= quants$class[2] ~ 'medium',
                              FC_worm > quants$class[2] ~ 'large'),
         .before = 'ID') %>% 
  write.xlsx('output/tables/worm_FC_PA_all.xlsx')

# save the full model 
worm_FC_PA_sig %>% 
  mutate(FC_class = case_when(FC_worm <= quants$class[1] ~ 'low',
                              FC_worm > quants$class[1] & 
                                FC_worm <= quants$class[2] ~ 'medium',
                              FC_worm > quants$class[2] ~ 'large'),
         .before = 'ID') %>% 
  write.xlsx('output/tables/worm_FC_PA_sig.xlsx')



worm_FC_PA_sig %>% 
  mutate(FC_class = case_when(FC_worm <= quants$class[1] ~ 'low',
                              FC_worm > quants$class[1] & FC_worm <= quants$class[2] ~ 'medium',
                              FC_worm > quants$class[2] ~ 'large'),
         .before = 'ID') %>% 
  ggplot(aes(x = fct_reorder(genome, FC_worm),
             y = FC_worm,
             color = FC_class)) +
  geom_point() +
  theme_classic() + 
  theme(
    axis.text.x = element_blank()
  )


worm_FC_PA_sig %>% 
  mutate(FC_class = ntile(FC_worm, 4),
         .before = 'ID')


# previous data
x_data_TOTAL <- read_csv("data/x_data_TOTAL.csv")

y_data_TOTAL <- read_csv("data/y_data_TOTAL.csv", 
                         col_names = FALSE)


x_data_TOTAL %>% 
  select(Broadphenotype, FC_class, phylogroup, Worm_metf_0) %>% 
  bind_cols(y_data_TOTAL) %>% 
  mutate(ID = paste0(FC_class,Worm_metf_0)) %>% 
  ggplot(aes(y = X1,
             x = fct_reorder(ID, X1),
             color = FC_class)) +
  geom_point() +
  theme(
    axis.text.x = element_blank()
  )


quants = x_data_TOTAL %>% 
  select(Broadphenotype, FC_class, phylogroup, Worm_metf_0) %>% 
  bind_cols(y_data_TOTAL) %>% 
  mutate(ID = paste0(FC_class,Worm_metf_0)) %>% 
  summarise(new_class = scales::percent(c(.20,.80)),
           class = quantile(X1, c(.20,.80)))

x_data_TOTAL %>% 
  select(Broadphenotype, FC_class, phylogroup, Worm_metf_0) %>% 
  bind_cols(y_data_TOTAL) %>% 
  mutate(ID = paste0(FC_class,Worm_metf_0)) %>% 
  mutate(new_class = case_when(X1 <= 1.35 ~ 'low',
                               X1 > 1.35 & X1 <= 2.19 ~ 'medium',
                               X1 > 2.19 ~ 'large')) %>% 
  ggplot(aes(y = X1,
             x = fct_reorder(ID, X1),
             color = new_class)) +
  geom_point() +
  theme(
    axis.text.x = element_blank()
  )







