
# libraries ---------------------------------------------------------------



library(tidyverse) # master library to deal with data frames
library(readxl) # read xlsx or xls files
library(FactoMineR) # for PCA
library(factoextra) # for PCA
library(here) # usefull to save plots in folders given a root
library(viridis) # color palette package
library(ComplexHeatmap) # yeah, complex heatmaps
library(openxlsx)
library(ggrepel)
library(patchwork)

theme_set(theme_light())

strain_db = read_excel("strain_db.xlsx")

# summarise phylogroups 
phylogr = strain_db %>%
  filter(Broadphenotype != 'Evolutionexperiment') %>%
  drop_na(Assembly) %>% 
  mutate(
    phylo_corrected = ifelse(phylogroup == 'Unknown', mash_group, phylogroup),
    phylo_corrected = factor(phylo_corrected, levels = c('A', 'B2', 'B1', 'C', 'D','E', 'F','G',
                                                    'cladeI', 'cladeII', 'cladeIII', 'cladeIV','cladeV',
                                                    'E or cladeI', 'fergusonii', 'albertii', 'Non Escherichia'
                                                    ))) %>%
  mutate(phylo_corrected = fct_explicit_na(phylo_corrected)) %>%
  group_by(phylo_corrected) %>%
  summarise(N = n())

# plot phylogroups
phylogr %>%
  mutate(Total = sum(N),
         Fraction = round((N/Total)*100,1),
         y_label_pos = Fraction + 3) %>%
  ggplot(aes(x = phylo_corrected, y = Fraction, fill = phylo_corrected)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(y = y_label_pos, label = Fraction), vjust = 1.6, size = 3.5) +
  labs(x = 'phylogroups', y = '% of total') +
  theme_light() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1))


dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'phylogrous_dist_barplot.pdf'),
             height = 8, width = 8, useDingbats = FALSE)



phylogr = strain_db %>%
  filter(Broadphenotype != 'Evolutionexperiment') %>%
  drop_na(Assembly) %>% 
  mutate(
    phylo_corrected = ifelse(phylogroup == 'Unknown', mash_group, phylogroup),
    phylo_corrected = factor(phylo_corrected, levels = c('A', 'B2', 'B1', 'C', 'D','E', 'F','G',
                                                         'cladeI', 'cladeII', 'cladeIII', 'cladeIV','cladeV',
                                                         'E or cladeI', 'fergusonii', 'albertii', 'Non Escherichia'
    ))) %>%
  mutate(phylo_corrected = fct_explicit_na(phylo_corrected)) %>%
  group_by(phylo_corrected, Broadphenotype) %>%
  summarise(N = n()) %>% ungroup

phylogr %>%
  mutate(Total = sum(N),
         Fraction = round((N/Total)*100,1),
         y_label_pos = Fraction + 3) %>%
  ggplot(aes(x = phylo_corrected, y = Fraction, fill = phylo_corrected)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(y = y_label_pos, label = Fraction), vjust = 1.6, size = 3.5) +
  labs(x = 'phylogroups', y = '% of total') +
  facet_wrap(~Broadphenotype) +
  theme_light() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1))

dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'phylogrous_phenotype_dist_barplot.pdf'),
             height = 10, width = 14, useDingbats = FALSE)




# antiSMASH results -------------------------------------------------------

antismash_summary = read_csv("antismashV5/antismash_results.csv") %>% 
  rename(Entry = X1,
         ID = Genome) %>% 
  separate(ID, c('ID', 'number'), sep = '_') %>% 
  select(-number)

library(readxl)
metadata = read_excel("D:/MRC_Postdoc/Pangenomic/Worm_imaging/MAIN_metadata.xlsx")

anti = antismash_summary %>% 
  mutate(product = str_replace(product, '\\[', ''),
         product = str_replace(product, '\\]', ''),
         product = str_replace_all(product, ' ', ''),
         product = str_replace_all(product, '\'', ''),
         kind = str_replace(kind, '\\[', ''),
         kind = str_replace(kind, '\\]', ''),
         kind = str_replace_all(kind, ' ', ''),
         kind = str_replace_all(kind, '\'', '')) %>% 
  separate_rows(product, sep = ',') %>% 
  full_join(metadata) %>% # full join to add the genomes without genome assembly
  drop_na(PG)


# how many genomes
anti %>% distinct(ID) %>% count
metadata %>% distinct(ID) %>% count


# how many products are encoded 
anti %>% 
  distinct(ID, product) %>% 
  count(product) %>% 
  drop_na(product) %>% 
  ggplot(aes(x = product, y = n)) +
  geom_bar(aes(fill =  product),stat='identity', color = 'black') +
  labs(x = 'Number of BGCs',
       y = 'Type of BGC',
       fill = 'BGC') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'antismash_products_summary.pdf'),
             height = 5, width = 8, useDingbats = FALSE)


# plot with percentages
anti %>% 
  distinct(ID, product) %>% 
  count(product) %>% 
  drop_na(product) %>%
  mutate(N = round(n/sum(414)*100,1)) %>% 
  ggplot(aes(x = product, y = N)) +
  geom_bar(aes(fill =  product),stat='identity', color = 'black') +
  labs(x = '% of BGCs',
       y = 'Type of BGC',
       fill = 'BGC') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'antismash_products_summary.pdf'),
             height = 6, width = 9, useDingbats = FALSE)









anti %>% filter(product == 'arylpolyene') %>% distinct(ID) %>% count


anti %>% 
  filter(ID == 'NT12010') %>%  View

# overwrite the main Metadata with the new info
# list of tables 
my_list = list('metadata' = metadata, 'metadata_antismash' = anti)


write.xlsx(my_list, "D:/MRC_Postdoc/Pangenomic/Worm_imaging/MAIN_metadata.xlsx", colNames = T, rowNames = T)




# pyseer WORM ------------------------------------------------------------------



# read pyseer data from worm phenotype


library(readr)
COGs_worm = read_delim("pyseer/COGs_worm_FC_lmm.tsv", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)








# manual manhattan plot
COGs_worm %>% 
  arrange(variant) %>% 
  mutate(log_lrt = -log10(`lrt-pvalue`),
         log_pval = -log10(`filter-pvalue`)) %>% 
  arrange(desc(af)) %>%
  mutate(variant_single = str_split(variant, pattern = '~~~') %>%  map_chr(., 1)) %>% 
  mutate(
        variant_2 = case_when(
                               grepl("group_",variant_single) ~ 'Unknown',
                               TRUE ~ variant_single),
         labels = case_when(log_lrt > 3 ~ variant_2),
        variant = factor(variant, levels = variant)) %>% 
  ggplot(aes(x = variant, y = log_lrt, color = af)) +
  geom_point() +
  labs(
    x = 'Genes',
    y = '-log10(lrt p-value)'
  ) +
  geom_hline(yintercept = c(0,3,6)) +
  geom_text_repel(aes(label = labels), box.padding = 1.1) +
  theme_classic() +
  theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 



dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'manhattan_pyseer_lrt_worm.pdf'),
             height = 10, width = 14, useDingbats = FALSE)







# manual manhattan plot
COGs_worm %>% 
  arrange(variant) %>% 
  mutate(log_lrt = -log10(`lrt-pvalue`),
         log_pval = -log10(`filter-pvalue`)) %>% 
  arrange(desc(af)) %>% 
  mutate(variant_single = str_split(variant, pattern = '~~~') %>%  map_chr(., 1)) %>% 
  mutate(variant_2 = case_when( variant_single == 'group_5967' ~ 'Group_5967',
                               grepl("group_",variant_single) ~ 'Unknown',
                               TRUE ~ variant_single),
         labels = case_when(log_pval > 6 ~ variant_2),
         variant = factor(variant, levels = variant)) %>% 
  ggplot(aes(x = variant, y = log_pval, color = af)) +
  geom_point() +
  labs(
    x = 'Genes',
    y = '-log10(p-value)'
  ) +
  geom_hline(yintercept = c(0,3,6)) +
  geom_text_repel(aes(label = labels), box.padding = 1.1) +
  theme_classic() +
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()) 


dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'manhattan_pyseer_pval_worm.pdf'),
             height = 10, width = 14, useDingbats = FALSE)








# pyseer BACT ------------------------------------------------------------------


# read pyseer data from worm phenotype


library(readr)
COGs_bact = read_delim("pyseer/COGs_bact_score_lmm.tsv", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)






# manual manhattan plot
COGs_bact %>% 
  arrange(desc(af)) %>% 
  mutate(log_lrt = -log10(`lrt-pvalue`),
         log_pval = -log10(`filter-pvalue`)) %>% 
  mutate(variant_single = str_split(variant, pattern = '~~~') %>%  map_chr(., 1)) %>% 
  mutate(variant_2 = case_when(grepl("group_",variant_single) ~ 'Unknown',
                               TRUE ~ variant_single),
         labels = case_when(log_lrt > 3 ~ variant_2),
         variant = factor(variant, levels = variant)) %>% 
  ggplot(aes(x = variant, y = log_lrt, color = af)) +
  geom_point() +
  labs(
    x = 'Genes',
    y = '-log10(lrt p-value)'
  ) +
  geom_hline(yintercept = c(0,3,6)) +
  scale_colour_gradient(low = "#4D1C00", high = '#FA5C00') +
  geom_text_repel(aes(label = labels), box.padding = 1.1) +
  theme_classic() +
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()) 



dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'manhattan_pyseer_lrt_bact.pdf'),
             height = 10, width = 14, useDingbats = FALSE)







# manual manhattan plot
COGs_bact %>% 
  arrange(desc(af)) %>% 
  mutate(log_lrt = -log10(`lrt-pvalue`),
         log_pval = -log10(`filter-pvalue`)) %>% 
  mutate(variant_single = str_split(variant, pattern = '~~~') %>%  map_chr(., 1)) %>% 
  mutate(variant_2 = case_when(grepl("group_",variant_single) ~ 'Unknown',
                               TRUE ~ variant_single),
         labels = case_when(log_pval > 3 ~ variant_2),
         variant = factor(variant, levels = variant)) %>% 
  ggplot(aes(x = variant, y = log_pval, color = af)) +
  geom_point() +
  labs(
    x = 'Genes',
    y = '-log10(p-value)'
  ) +
  geom_hline(yintercept = c(0,3,6)) +
  geom_text_repel(aes(label = labels), box.padding = 1.1) +
  scale_colour_gradient(low = "#4D1C00", high = '#FA5C00') +
  theme_classic() +
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()) 


dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'manhattan_pyseer_pval_bact.pdf'),
             height = 10, width = 14, useDingbats = FALSE)




# NEW pyseer worms --------------------------------------------------------------



# get the data

gene_presence_absence <- read_delim("D:/MRC_Postdoc/Pangenomic/phylo/original_data/panaroo_clean/gene_presence_absence.Rtab", 
                                    "\t", escape_double = FALSE, trim_ws = TRUE)

summary.wb = read_excel("Summary_worm_bacteria.xlsx") %>% select(-`...1`)


# create df
genes = COGs_worm$variant

summary_FC_bygene = data.frame(genes)
summary_FC_bygene['FC'] = 0

# # # # # # # # # # 
## RUN ONLY ONCE ##
# # # # # # # # # #
# Just a note, I removed all the strains that were not suitable 
# from the worm experiments as well

for (i in 1:length(genes)){
  
  print(paste('gene ', as.character(i), 'out of ', as.character(length(genes))))
  
  cosa = gene_presence_absence %>%
    filter(Gene == genes[i]) %>%
    pivot_longer(cols = !Gene, names_to = 'ID', values_to = 'presence') %>%
    separate(ID, into = c('ID', 'number'), sep = '_') %>% 
    select(-number) %>%
    left_join(summary.wb) %>%
    filter(Worm_metf_0 < 4000,
           biofilm_50mM == 'normal') %>%   # FILTERS
    mutate(presence = as.factor(presence))%>% 
    group_by(presence) %>% 
    summarise(mean = mean(FC_worm)) %>% 
    pivot_wider(names_from = presence, values_from = mean) %>% 
    mutate(FC = `1`/`0`) %>% 
    select(FC)
  
  summary_FC_bygene[summary_FC_bygene$gene == genes[i],]$FC = cosa$FC
  
}

summary_FC_bygene = summary_FC_bygene %>% tibble()


genes_plot = c('cra_2', 'scrY', 'scrB', 'scrK')
# add this to filter genes "& variant_2 %in% genes_plot"
# generate dataset 
COGs_plot = COGs_worm %>% 
  arrange(variant) %>% 
  left_join(summary_FC_bygene %>% 
              select(variant = genes,
                     FC)) %>% 
  mutate(log_lrt = -log10(`lrt-pvalue`),
         log_pval = -log10(`filter-pvalue`)) %>% 
  arrange(desc(af)) %>%
  mutate(variant_single = str_split(variant, pattern = '~~~') %>%  map_chr(., 1)) %>% 
  mutate(
    variant_2 = case_when(
      grepl("group_",variant_single) ~ 'Unknown',
      TRUE ~ variant_single),
    point_shape = case_when(variant_2 == 'Unknown' ~ 2,
                            TRUE ~ 1),
    point_shape = as.factor(point_shape),
    labels = case_when(log_pval > 5 & variant_2 != 'Unknown'  ~ variant_2),
    variant = factor(variant, levels = variant),
    FC_class = cut(FC, breaks=c(-10,1,10), labels=c("negative", "positive")),
    FC_class = as.factor(FC_class))



# line plot
p1 = COGs_plot %>%
  ggplot(aes(x = variant, y = af, color = af, group = 1)) +
  geom_line(size = 1) +
  labs(x = 'Genes',
       y = 'Gene\nfrequency') + 
  theme_classic() +
  theme( axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         legend.position = "none",
         axis.title.y = element_text(size = 16),
         axis.title.x = element_text(size = 16))

# manhattan plot
p2 = COGs_plot %>% 
  ggplot(aes(x = variant, y = log_pval, color = FC_class)) +
  geom_point(aes(shape = point_shape)) +
  labs(
    x = '',
    y = '-log10(p-value)'
  ) +
  ylim(0,10.1) +
  geom_hline(yintercept = c(0,5)) +
  geom_text_repel(aes(label = labels),
                  nudge_y = 2,
                  # direction    = "x",
                  box.padding = 0.3, 
                  vjust        = 0,
                  segment.size = 0.2,
                  segment.color = "grey90",
                  segment.alpha = 0.5,
                  size = 5)  +
  theme_classic() +
  scale_color_manual(values = c("#00AFBB", 
                                "#FC4E07"),
                     name = 'Effect size\ndirection') +
  scale_shape_discrete(name = 'Gene type', labels = c('Known', 'Unknown')) +
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(), 
    axis.title.y = element_text(size = 16),
    legend.position = "right",
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 13)) 

# pathworks
p2 / p1 +
plot_layout(ncol = 1, heights = c(13, 1))

dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'manhattan_pyseer_pval_worm_CECAD.pdf'),
             height = 8, width = 11, useDingbats = FALSE)





# NEW pyseer bac ----------------------------------------------------------




# generate dataset 
COGs_plot = COGs_bact %>% 
  arrange(variant) %>% 
  left_join(summary_FC_bygene %>% 
              select(variant = genes,
                     FC)) %>% 
  mutate(log_lrt = -log10(`lrt-pvalue`),
         log_pval = -log10(`filter-pvalue`)) %>% 
  arrange(desc(af)) %>%
  mutate(variant_single = str_split(variant, pattern = '~~~') %>%  map_chr(., 1)) %>% 
  mutate(
    variant_2 = case_when(
      grepl("group_",variant_single) ~ 'Unknown',
      TRUE ~ variant_single),
    point_shape = case_when(variant_2 == 'Unknown' ~ 2,
                            TRUE ~ 1),
    point_shape = as.factor(point_shape),
    labels = case_when(log_pval > 5 & variant_2 != 'Unknown' ~ variant_2),
    variant = factor(variant, levels = variant),
    FC_class = cut(FC, breaks=c(-10,1,10), labels=c("negative", "positive")),
    FC_class = as.factor(FC_class))



# line plot
p1 = COGs_plot %>%
  ggplot(aes(x = variant, y = af, color = af, group = 1)) +
  geom_line(size = 1) +
  labs(x = 'Genes',
       y = 'Gene\nfrequency') + 
  theme_classic() +
  theme( axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         legend.position = "none",
         axis.title.y = element_text(size = 16),
         axis.title.x = element_text(size = 16))

# manhattan plot
p2 = COGs_plot %>% 
  ggplot(aes(x = variant, y = log_pval, color = FC_class)) +
  geom_point(aes(shape = point_shape)) +
  labs(
    x = '',
    y = '-log10(p-value)'
  ) +
  geom_hline(yintercept = c(0,5)) +
  ylim(0,10.1) +
  geom_text_repel(aes(label = labels),
                  nudge_y = 2,
                  # direction    = "x",
                  box.padding = 0.1, 
                  vjust        = 0,
                  segment.size = 0.2,
                  segment.color = "grey90",
                  segment.alpha = 0.5,
                  size = 3)  +
  theme_classic() +
  scale_color_manual(values = c("#00AFBB", 
                                "#FC4E07"),
                     name = 'Effect size\ndirection') +
  scale_shape_discrete(name = 'Gene type', labels = c('Known', 'Unknown')) +
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(), 
    axis.title.y = element_text(size = 16),
    legend.position = "right",
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 13)) 

# pathworks
p2 / p1 +
  plot_layout(ncol = 1, heights = c(13, 1))

dev.copy2pdf(device = cairo_pdf,
             file = here('R_plots', 'manhattan_pyseer_pval_bact_final.pdf'),
             height = 8, width = 11, useDingbats = FALSE)






# ANN data preparation ----------------------------------------------------

# here I'll prepare the data to see if I can fit a simple neural network

# first let's work the gene matrix
# create a copy just in case
gene_presence_absence_alt = gene_presence_absence

# remove singletons
gene_presence_absence_alt = gene_presence_absence_alt %>% 
  rowwise(Gene) %>% 
  mutate( suma = sum(c_across(where(is.numeric))), .before = NT12001_189) %>% 
  filter(suma > 1) %>% 
  select(-suma)

# make the transpose
gene_presence_absence_alt = gene_presence_absence_alt %>%
  gather(key = ID, value = val, 2:ncol(gene_presence_absence_alt)) %>% 
  spread(Gene, val)

gene_presence_absence_alt = gene_presence_absence_alt %>% 
  mutate(ID = str_split(ID, pattern = '_', simplify = TRUE)[,1])





neural_net_test = summary.wb %>% 
  filter(biofilm_50mM == 'normal',
         Worm_metf_0 < 4000) %>% 
  drop_na(phylogroup) %>% 
  arrange(ID) %>% 
  distinct(ID, .keep_all = T) %>% 
  select(-(SD_Bact_metf_0:SD_Bact_metf_200),-Strainname,-PG,-Well,-B_phenotype,
         -biofilm_0mM,-biofilm_50mM,-log2FC_worm) %>% 
  replace_na(list(Broadphenotype = 'Unknown')) %>% 
  left_join(gene_presence_absence_alt) %>% 
  select(-ID)


y_data = neural_net_test$FC_worm %>% write.table(here('NN','y_data.csv'),quote=F, col.names = F, row.names = F)
x_data = neural_net_test %>% select(-FC_worm) %>% write_csv(here('NN','x_data.csv'))



