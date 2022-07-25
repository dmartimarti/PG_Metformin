

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
  left_join(gene_pa_all)

# case for only genes that are significant
worm_FC_PA_sig = metadata %>% 
  mutate(genome = str_sub(fasta, 1, -7),
         .before= 'ID') %>% 
  filter(genome %in% input_pyseer$IDs) %>% 
  distinct(genome, .keep_all = TRUE) %>% 
  filter(Discard == 'No') %>% 
  left_join(input_pyseer %>% 
              rename(genome = IDs)) %>% 
  left_join(gene_pa_sig)
# save the full model 
worm_FC_PA_all %>% 
  write.xlsx('output/tables/worm_FC_PA_all.xlsx')

# save the full model 
worm_FC_PA_sig %>% 
  write.xlsx('output/tables/worm_FC_PA_sig.xlsx')


