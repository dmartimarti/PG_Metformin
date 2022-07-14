# libraries

library(tidyverse)
library(readr)
library(here)
library(cowplot)
library(rstatix)
library(ggrepel)





# folders and data loading ------------------------------------------------

folders = c('rcdA_MvsC', 'rcdA_M_vs_OP50_M', 'rcdA_C_vs_OP50_C', 
            'OP50_MvsC', 'DM_MvsC', 'DM_C_vs_OP50_C', 'DM_C_vs_CRP_C',
            'CRP_MvsC')

OP50_MvC = read_delim("rcdA_MvsC/enrichment.tsv", 
                                   "\t", escape_double = FALSE, trim_ws = TRUE, 
                                   skip = 1) %>% 
  rename(test_type = `Categorical column tested against the 'selection'`,
         direction = `Value in the 'selection' column`,
         term = `Term from this categorical column.`,
         p_val = `P value resulting from the Fisher exact test.`,
         fdr = `Benjamini Hochberg false discovery rate.`)




tests = unique(OP50_MvC$test_type)

# I can loop here

for (i in 1:length(tests)) {

  p = OP50_MvC %>% 
    filter(test_type == tests[1]) %>% 
    mutate(category = cut(fdr, 
                          breaks=c(-Inf, 0.001,  0.01, 0.05, Inf), 
                          labels=c("< 0.001","< 0.01","<0.05", 'NS') )) %>% 
    ggplot(aes(y = term, x = direction, fill = category)) +
    geom_tile() +
    labs(title = tests[1]) +
    scale_fill_manual(values = c('#2E2ED1', '#6B6BD1', '#A3A3D1', 'white')) +
    theme_classic()
  
  ggsave(plot = p, filename = here('rcdA_MvsC', paste0('enrich_',tests[1],'.pdf')), 
         device = 'pdf', height = 8, width = 10)
}








# automatization of everything --------------------------------------------



# only the folders that contain the file
# ideally I should do a try/except for the loop, but it works fine for now
folders = c('rcdA_MvsC', 'rcdA_M_vs_OP50_M', 'OP50_MvsC', 'DM_C_vs_OP50_C', 'CRP_MvsC')

for (folder in folders){
  print(paste0('printing plots from folder ', folder))
  enrich = read_delim(paste0(folder,"/enrichment.tsv"), 
                        "\t", escape_double = FALSE, trim_ws = TRUE, 
                        skip = 1) %>% 
    rename(test_type = `Categorical column tested against the 'selection'`,
           direction = `Value in the 'selection' column`,
           term = `Term from this categorical column.`,
           p_val = `P value resulting from the Fisher exact test.`,
           fdr = `Benjamini Hochberg false discovery rate.`)

  tests = unique(enrich$test_type)
  
  # I can loop here
  
  for (i in 1:length(tests)) {
    
    p = enrich %>% 
      filter(test_type == tests[i]) %>% 
      mutate(`P-value` = cut(fdr, 
                            breaks=c(-Inf, 0.001,  0.01, 0.05, Inf), 
                            labels=c("< 0.001","< 0.01","<0.05", 'NS') )) %>% 
      ggplot(aes(y = term, x = direction, fill = `P-value`)) +
      geom_tile(width = 1, height = 1) +
      labs(title = tests[i],
           x = 'Direction',
           y = 'Ontology term') +
      scale_fill_manual(values = c('#2E2ED1', '#6B6BD1', '#A3A3D1', 'white')) +
      # scale_fill_manual(values = c('grey30', 'grey50', 'grey70', 'white')) +
      theme_cowplot(15)
    
    ggsave(plot = p, filename = here(folder, paste0('enrich_',tests[i],'.pdf')), 
           device = 'pdf', height = 8, width = 7)
  }
}





# transporters from CRP and OP50 ------------------------------------------




  
  
## files from Perseus -----------


op50_diff <- read_delim("prot_data_perseus/op50_MvC_prot_diff.tsv", 
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)

op50_diff = op50_diff %>% 
  pivot_longer(cols = OP50_0mM_C1:OP50_50mM_C5, 
               names_to = 'sample', values_to = 'intensity') %>% 
  select(sample, intensity, 
         gene = `Gene names`,
         protein = `Protein names`,
         kegg = `KEGG brite name`,
         everything(),
         -`Razor + unique peptides`,
         -Intensity,
         -`Majority protein IDs`,
         -`MS/MS count`) %>% 
  separate(sample, into = c('sample','metformin','replicate'),
           sep = '_') %>% 
  mutate(Cluster = case_when(Cluster == 'Cluster -1326' ~ 'Up',
                             Cluster == 'Cluster -1325' ~ 'Down'))




crp_diff <- read_delim("prot_data_perseus/crp_MvsC_prot_diff.tsv", 
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)

crp_diff = crp_diff %>% 
  pivot_longer(cols = CRP_0mM_C1:CRP_50mM_C5, 
               names_to = 'sample', values_to = 'intensity') %>% 
  select(sample, intensity, 
         gene = `Gene names`,
         protein = `Protein names`,
         kegg = `KEGG brite name`,
         everything(),
         -`Razor + unique peptides`,
         -Intensity,
         -`Majority protein IDs`,
         -`MS/MS count`) %>% 
  separate(sample, into = c('sample','metformin','replicate'),
           sep = '_') %>% 
  mutate(Cluster = case_when(Cluster == 'Cluster -634' ~ 'Up',
                             Cluster == 'Cluster -633' ~ 'Down'))



op50_diff %>% 
  separate_rows(kegg, sep = ';') %>% 
  filter(kegg == 'Transporters') %>% 
  drop_na(Cluster) %>% 
  distinct(gene,.keep_all = T) %>% 
  select(-replicate,-intensity,-id,-sample,-metformin) %>% 
  select(gene, protein, direction = Cluster, everything()) %>% 
  arrange(direction) %>% 
  write_csv('transporters_OP50.csv')

crp_diff %>% 
  separate_rows(kegg, sep = ';') %>% 
  filter(kegg == 'Transporters') %>% 
  drop_na(Cluster) %>% 
  distinct(gene,.keep_all = T) %>% 
  select(-replicate,-intensity,-id,-sample,-metformin) %>% 
  select(gene, protein, direction = Cluster, everything()) %>% 
  arrange(direction) %>% 
  write_csv('transporters_CRP.csv')


# genes: mdtK and ybjK

genes_of_interest = c('mdtK', 'ybjK','ybjJ')

op50_diff %>% 
  filter(gene %in% genes_of_interest) %>% 
  ggplot(aes(x = metformin, y = intensity, fill = metformin)) +
  geom_boxplot(show.legend = F) +
  geom_point(position = position_jitterdodge()) +
  cowplot::theme_cowplot(14) +
  facet_wrap(~gene)

ggsave('mdtK_OP50.pdf', 
       height = 5, width = 6)



crp_diff %>% 
  filter(gene %in% genes_of_interest) %>% 
  ggplot(aes(x = metformin, y = intensity, fill = metformin)) +
  geom_boxplot(show.legend = F) +
  geom_point(position = position_jitterdodge()) +
  cowplot::theme_cowplot(14) +
  facet_wrap(~gene)

ggsave('mdtK_CRP.pdf', 
       height = 5, width = 6)




# nice volcano plots ####

## OP50 ---------


op50 = read_delim("prot_data_perseus/op50_MvC_prot.tsv", 
                  delim = "\t", escape_double = FALSE, 
                  trim_ws = TRUE)


op50 = op50 %>% 
  pivot_longer(cols = OP50_0mM_C1:OP50_50mM_C5, 
               names_to = 'sample', values_to = 'intensity') %>% 
  select(sample, intensity, 
         gene = `Gene names`,
         protein = `Protein names`,
         kegg = `KEGG brite name`,
         everything(),
         -`Razor + unique peptides`,
         -Intensity,
         -`Majority protein IDs`,
         -`MS/MS count`) %>% 
  separate(sample, into = c('sample','metformin','replicate'),
           sep = '_')

op50_stats = op50 %>% 
  mutate(metformin = as.factor(metformin)) %>% 
  group_by(sample, gene, protein, kegg) %>% 
  t_test(intensity ~ metformin, p.adjust.method = 'fdr', detailed = T) 


# volcano plot
pval = -log10(0.05)
estimate_thr = 1

op50_stats %>% 
  drop_na(p) %>% 
  mutate(estimate = estimate*-1) %>% # reorder log2FC for treat vs control
  mutate(logP = -log10(p)) %>% 
  mutate(groups = case_when(logP > pval & estimate > estimate_thr ~ 'up_reg',
                            logP > pval & estimate < -estimate_thr ~ 'down_reg',
                            TRUE ~ 'non_sig'),
         gene_plot = case_when(logP > pval & abs(estimate) > estimate_thr ~ gene)) %>% 
  ggplot(aes(x = estimate, y = logP)) +
  geom_point(aes(color = groups), show.legend = F) +
  geom_text_repel(aes(label = gene_plot,box.padding = 0.1,
                      max.overlaps = 10)) +
  scale_color_manual(values = c('blue', 'gray90', 'red')) +
  labs(
    x = 'Log2FC',
    y = '-log10(P-value)'
  ) +
  theme_cowplot(15)

ggsave('VOLCANO_PLOTS/OP50_volcano.pdf',
       height = 7, width = 7)

# test gnes to 
op50 %>% 
  filter(gene == 'yiaY') %>% 
  ggplot(aes(x = metformin, y = intensity)) +
  geom_boxplot()


## OP50 ---------


crp = read_delim("prot_data_perseus/crp_MvsC_prot.tsv", 
                  delim = "\t", escape_double = FALSE, 
                  trim_ws = TRUE)


crp = crp %>% 
  pivot_longer(cols = CRP_0mM_C1:CRP_50mM_C5, 
               names_to = 'sample', values_to = 'intensity') %>% 
  select(sample, intensity, 
         gene = `Gene names`,
         protein = `Protein names`,
         kegg = `KEGG brite name`,
         everything(),
         -`Razor + unique peptides`,
         -Intensity,
         -`Majority protein IDs`,
         -`MS/MS count`) %>% 
  separate(sample, into = c('sample','metformin','replicate'),
           sep = '_')

crp_stats = crp %>% 
  mutate(metformin = as.factor(metformin)) %>% 
  group_by(sample, gene, protein, kegg) %>% 
  t_test(intensity ~ metformin, p.adjust.method = 'fdr', detailed = T) 


# volcano plot
pval = -log10(0.05)
estimate_thr = 1

crp_stats %>% 
  drop_na(p) %>% 
  mutate(estimate = estimate*-1) %>% # reorder log2FC for treat vs control
  mutate(logP = -log10(p)) %>% 
  mutate(groups = case_when(logP > pval & estimate > estimate_thr ~ 'up_reg',
                            logP > pval & estimate < -estimate_thr ~ 'down_reg',
                            TRUE ~ 'non_sig'),
         gene_plot = case_when(logP > pval & abs(estimate) > estimate_thr ~ gene)) %>% 
  ggplot(aes(x = estimate, y = logP)) +
  geom_point(aes(color = groups), show.legend = F) +
  geom_text_repel(aes(label = gene_plot), box.padding = 0.1,
                  max.overlaps = 10) +
  scale_color_manual(values = c('blue', 'gray90', 'red')) +
  labs(
    x = 'Log2FC',
    y = '-log10(P-value)'
  ) +
  theme_cowplot(15)

ggsave('VOLCANO_PLOTS/CRP_volcano.pdf',
       height = 7, width = 7)

# test gnes to 
crp %>% 
  filter(gene == 'asr') %>% 
  ggplot(aes(x = metformin, y = intensity)) +
  geom_boxplot()



# plot with REAL volcano!!

require(magick)

volcano_plot = crp_stats %>% 
  drop_na(p) %>% 
  mutate(estimate = estimate*-1) %>% # reorder log2FC for treat vs control
  mutate(logP = -log10(p)) %>% 
  mutate(groups = case_when(logP > pval & estimate > estimate_thr ~ 'up_reg',
                            logP > pval & estimate < -estimate_thr ~ 'down_reg',
                            TRUE ~ 'non_sig'),
         gene_plot = case_when(logP > pval & abs(estimate) > estimate_thr ~ gene)) %>% 
  ggplot(aes(x = estimate, y = logP)) +
  geom_point(aes(color = groups), show.legend = F) +
  geom_text_repel(aes(label = gene_plot), box.padding = 0.1,
                  max.overlaps = 10) +
  scale_color_manual(values = c('blue', 'gray90', 'red')) +
  labs(
    x = 'Log2FC',
    y = '-log10(P-value)'
  ) +
  theme_cowplot(15)


ggdraw() +
  draw_image("https://easydrawingguides.com/wp-content/uploads/2018/09/Volcano-10.png",
             y = -0.5, x = 0.02) + 
  draw_plot(volcano_plot) +
  theme(plot.margin = margin(0,0,2.5,0, "cm"))





