# libraries --------------------------------------------------------------

library(tidyverse)
library(readr)
library(readxl)
theme_set(theme_light())
library(openxlsx)
library(here)
library(ggrepel)
library(viridis)
library(openxlsx)
# library(tidymodels)
library(cowplot)

# load data ---------------------------------------------------------------

keio.rat = read_excel("D:/MRC_Postdoc/Pangenomic/biolog/KEIO/fixed_data/KEIO_output/summary/AUC_KEIO_screen.xlsx", 
                      sheet = "ratios") %>% select(-`...1`)

aska.rat = read_excel("D:/MRC_Postdoc/Pangenomic/biolog/ASKA/fixed_data/Output/summary/AUC_raw_ASKAscreen.xlsx", 
                      sheet = "ratios_by_col") %>% select(-`...1`)

# select data
keio.rat = keio.rat %>% select(Strain, keio_mean = Mean_correct)
aska.rat = aska.rat %>% select(Strain, aska_mean = Mean_correct)

merg_data = keio.rat %>% left_join(aska.rat) %>% drop_na(keio_mean, aska_mean)


# data transformation -----------------------------------------------------
# the idea for this part is to calculate the center of the distributions and 
# then calculate the distance of each point from the center. 
# to add labels that might get lost, I also add the labels for the top and bottom
# 40 hits (number can change) from both distributions

# mean as the center of the circle
x.1 = keio.rat %>% summarise(Mean = mean(keio_mean))
y.1 = aska.rat %>% summarise(Mean = mean(aska_mean))

# data transformations to separate hits
merg_data = merg_data %>% 
  mutate(x = x.1$Mean,
         y = y.1$Mean,
         distance = sqrt((x - keio_mean)**2 + (y - aska_mean)**2),
         keio_2 = keio_mean - x,
         aska_2 = aska_mean - y,
         distance2 = sqrt(keio_2**2 + aska_2**2),
         ratios = keio_mean/aska_mean,
         lograt = log10(ratios),
         rest = abs(keio_mean - aska_mean)) 

# get top resistant/sensitive hits
n = 40
bot_keio = merg_data %>% arrange(keio_mean) %>% head(n) %>% select(Strain) %>% t %>% as.character()
top_keio = merg_data %>% arrange(desc(keio_mean)) %>% head(n) %>% select(Strain) %>% t %>% as.character()
bot_aska = merg_data %>% arrange(aska_mean) %>% head(n) %>% select(Strain) %>% t %>% as.character()
top_aska = merg_data %>% arrange(desc(aska_mean)) %>% head(n) %>% select(Strain) %>% t %>% as.character()

lab = unique(c(bot_aska, bot_keio, top_aska, top_keio))

# label the 
merg_data = merg_data %>% 
  mutate(labels = case_when(rest > 0.38 ~ Strain,
                            Strain %in% lab ~ Strain))


merg_data %>% 
  # filter(!(Strain %in% c('yqaA', 'yajC'))) %>% 
  ggplot(aes(x = aska_mean, y = keio_mean, size = rest, colour = distance)) +
  geom_hline(yintercept = 0.42, color = 'black', size = 1.2, alpha = 0.4) +
  geom_vline(xintercept = 0.49, color = 'black', size = 1.2, alpha = 0.4) +
  geom_point() +
  geom_text_repel(aes(label = labels),
                  max.overlaps = Inf) +
  scale_color_viridis_c(name = 'Distance', option = 'viridis', direction = -1) +
  labs(x = 'ASKA mean ratio',
       y = 'KEIO mean ratio',
       size = "Ratio \ndifference") +
  theme_half_open(17)


merg_data %>% 
  # filter(!(Strain %in% c('yqaA', 'yajC'))) %>% 
  ggplot(aes(x = aska_2, y = keio_2, size = rest, colour = distance2)) +
  geom_hline(yintercept = 0, color = 'black', size = 1.2, alpha = 0.4) +
  geom_vline(xintercept = 0, color = 'black', size = 1.2, alpha = 0.4) +
  geom_point() +
  geom_text_repel(aes(label = labels),
                  max.overlaps = Inf) +
  scale_color_viridis_c(name = 'Distance', option = 'viridis', direction = -1) +
  labs(x = 'ASKA mean ratio',
       y = 'KEIO mean ratio',
       size = "Ratio \ndifference") +
  # guides(size = 'none') +
  theme_cowplot(17)

ggsave(file = here('summary','KEIO_ASKA_summary.pdf'),
       width = 13, height = 10)  


ggsave(file = here('summary','KEIO_ASKA_summary_v2.pdf'),
       width = 10, height = 8)  








# Tukey test --------------------------------------------------------------
library(ODWGtools)

# I am selecting the top and bottom genes from each screen by calculating the Tukey outliers
# from each distribution

# get outliers from keio
keio.out = keio.rat %>% 
  mutate(outliers = tukey_outliers(.$keio_mean)) %>% 
  filter(outliers != 'not outlier') %>% 
  mutate(direction = case_when(keio_mean > 0.5 ~ 'top',
                              keio_mean <= 0.5 ~ 'bottom')) 
# get outliers from aska
aska.out = aska.rat %>% 
  mutate(outliers = tukey_outliers(.$aska_mean)) %>% 
  filter(outliers != 'not outlier') %>% 
  mutate(direction = case_when(aska_mean > 0.5 ~ 'top',
                               aska_mean <= 0.5 ~ 'bottom')) 

# how many outliers in each part
keio.out %>% count(direction)

aska.out %>% count(direction)


keio.out %>% filter(direction == 'bottom') %>% 
  select(Strain) %>% t %>% as.character%>% 
  write.table('bottom_keio.txt', quote = F, row.names = F, col.names = F)

keio.out %>% filter(direction == 'top') %>% 
  select(Strain) %>% t %>% as.character%>% 
  write.table('top_keio.txt', quote = F, row.names = F, col.names = F)

aska.out %>% filter(direction == 'bottom') %>% 
  select(Strain) %>% t %>% as.character%>% 
  write.table('bottom_aska.txt', quote = F, row.names = F, col.names = F)

aska.out %>% filter(direction == 'top') %>% 
  select(Strain) %>% t %>% as.character%>% 
  write.table('top_aska.txt', quote = F, row.names = F, col.names = F)



# alternative way of top/bottom genes, using a threshold

# get outliers from keio
# keio.out_b = keio.rat %>%
#   arrange(keio_mean) %>% 
#   head(100) %>% 
#   mutate(direction = 'bottom')
# 
# keio.out_t = keio.rat %>% 
#   arrange(desc(keio_mean)) %>% 
#   head(100) %>% 
#   mutate(direction = 'top')
# 
# keio.out = rbind(keio.out_b, keio.out_t)
# 
# aska.out_b = aska.rat %>%
#   arrange(aska_mean) %>% 
#   head(100) %>% 
#   mutate(direction = 'bottom')
# 
# aska.out_t = aska.rat %>% 
#   arrange(desc(aska_mean)) %>% 
#   head(100) %>% 
#   mutate(direction = 'top')
# 
# aska.out = rbind(aska.out_b, aska.out_t)
# 
# 


# enrichment  -------------------------------------------------------------

### GO enrichment

# load ecocyc raw data
ecocyc = read_delim("All_genes_+_pathways_+_GO,_E.coli_K12_common_names.txt", 
                                                            "\t", escape_double = FALSE, trim_ws = TRUE)
# e coli proteome, a bypass for gene name synonyms
ecoprot = read_delim("uniprot-proteome_UP000000625.txt", 
                     "\t", escape_double = FALSE, trim_ws = TRUE)


# column names are bit weird
names(ecocyc)


ecocyc = ecocyc %>% 
  rename(gene = `Gene Name`,
         pathway = `Pathways of gene`,
         promDown = `Promoters - promoters regulated by gene`,
         promUP = `Promoters - promoters upstream of gene`,
         GO_BP = `GO terms (biological process)`,
         GO_CC = `GO terms (cellular component)`,
         GO_MF = `GO terms (molecular function)`)


ecoprot = ecoprot %>% 
  rename(gene_names = `Gene names`)


# super ugly and SLOW loop to get synonyms from the E. coli proteome
ecocyc['gene_names'] = ''

for (i in 1:dim(ecocyc)[1]){
  for (j in 1:dim(ecoprot)){
    if (str_detect(ecoprot$gene_names[j], ecocyc$gene[i]) == TRUE){
      ecocyc$gene_names[i] = ecoprot$gene_names[j]
    }
  } 
  print(i/dim(ecocyc)[1] * 100)
}

ecocyc = ecocyc %>% 
  mutate(gene_names = case_when(gene_names == "" ~ gene,
                                gene_names != "" ~ gene_names))

eco_path = ecocyc %>% 
  select(gene, gene_names, pathway) %>% drop_na(pathway) %>% 
  separate_rows(pathway, sep = " // ") %>% 
  separate_rows(gene_names, sep = ' ')
  

eco_path = eco_path %>% select(gene = gene_names, pathway)

# subset the complete database for the genes you have in each dataset
eco_path.keio = eco_path %>% filter(gene %in% keio.rat$Strain)
eco_path.aska = eco_path %>% filter(gene %in% aska.rat$Strain)


# save eco_path.aska as the main table for next analyses with enrichment
write_csv(eco_path.aska, here('summary', 'EcoCyc_genes_pathways.csv'))



# Hypergeometric function -------------------------------------------------


# hypergeometric test function

enrich = function(gene, db){
  # initiate variables
  pval = c()
  m_total = c()
  x_total = c()
  k_total = c()
  gene_in_cat = c()
  db = as.data.frame(db)
  cats = unique(db[,2])
  
  for (cat in cats){
    subcat = db[db[,2] == cat,]
    N = (db %>% distinct(.[,1]) %>% count())$n
    m = dim(subcat)[1]
    n = N - m
    x = sum(gene %in% subcat[,1])
    k = sum(gene %in% db[,1]) # genes with at least 1 annotation!
    p = phyper(q=x-1, m=m, n=n, k=k, lower.tail=FALSE)
    
    # save variables
    m_total = c(m_total, m)
    x_total = c(x_total, x)
    k_total = c(k_total, k)
    gene_in_cat = c(gene_in_cat)
    pval = c(pval, p)
  }
  
  # build the table
  table = tibble(categories = cats, N = N, elm_in_cat = m_total, gene_in_cat = x_total, k_tot = k, pval = pval) %>% 
    mutate(
      p.stars = gtools::stars.pval(pval),
      fdr = p.adjust(pval, method = 'fdr'),
      fdr.stars = gtools::stars.pval(fdr)
           ) %>% 
    arrange(pval)
  
  return(table)
  
}

# calculate enrichment for the different sets
# keio top
genes = (keio.out %>% filter(direction == 'top'))$Strain
keio.top.enrich = enrich(genes, db = eco_path.keio) %>% mutate(direction = 'top')

# keio bottom
genes = (keio.out %>% filter(direction == 'bottom'))$Strain
keio.bot.enrich = enrich(genes, db = eco_path.keio) %>% mutate(direction = 'bottom')

# aska top
genes = (aska.out %>% filter(direction == 'top'))$Strain
aska.top.enrich = enrich(genes, db = eco_path.aska) %>% mutate(direction = 'top')

# keio bottom
genes = (aska.out %>% filter(direction == 'bottom'))$Strain
aska.bot.enrich = enrich(genes, db = eco_path.aska) %>% mutate(direction = 'bottom')

# save lists in an excel file
# save info


list_of_datasets = list('Keio EcoCyc' = keio.top.enrich %>% bind_rows(keio.bot.enrich) %>% filter(pval <= 0.1), 
                        'Aska EcoCyc' = aska.top.enrich %>% bind_rows(aska.bot.enrich) %>% filter(pval <= 0.1))

write.xlsx(list_of_datasets, here('summary','Ecocyc_enrichment.xlsx'), colNames = T, rowNames = T) 

# plot
expanded = keio.top.enrich %>% filter(pval <= 0.05) %>% mutate(direction = 'top') %>% 
  bind_rows(keio.bot.enrich %>% filter(pval <= 0.05) %>% mutate(direction = 'bottom')) %>% 
  select(categories, pval, p.stars, direction) %>% 
  expand(categories, direction)


expanded %>% left_join(keio.top.enrich %>% filter(pval <= 0.05) %>% mutate(direction = 'top') %>% 
                         bind_rows(keio.bot.enrich %>% filter(pval <= 0.05) %>% mutate(direction = 'bottom')) %>% 
                         select(categories, pval, p.stars, direction)) %>% 
  replace_na(list(pval = 1, p.stars = ''))



# plots -------------------------------------------------------------------


# plot KEIO
list_of_datasets$`Keio EcoCyc` %>% 
  filter(fdr  <= 0.1) %>% 
  mutate(Category = cut(fdr , 
                        breaks=c(-Inf, 0.001,  0.01, 0.05, 0.1, Inf), 
                        labels=c("<0.001","<0.01","<0.05", '<0.1', 'NS') )) %>% 
  mutate(direction = case_when(direction == 'bottom' ~ 'Sensitive',
                               direction == 'top' ~ 'Resistant')) %>% 
  ggplot(aes(y = categories, x = direction, fill = Category)) +
  geom_tile() +
  labs(
    x = 'Resistance',
    y = 'Terms'
  ) +
  scale_fill_manual(values = c('#2229F0', '#2292F0', '#22E4F0', '#C7F0F0', 'white'),
                    name = 'P-val') +
  theme_classic()

ggsave(here('summary', 'keio_EcoCyc_enrichment_pval.pdf'), height =  7, width = 8)



# plot ASKA
list_of_datasets$`Aska EcoCyc` %>% 
  filter(fdr <= 0.05) %>% 
  mutate(Category = cut(fdr, 
                        breaks=c(-Inf, 0.001,  0.01, 0.05, 0.1, Inf), 
                        labels=c("<0.001","<0.01","<0.05", '<0.1', 'NS') )) %>% 
  mutate(direction = case_when(direction == 'bottom' ~ 'Sensitive',
                               direction == 'top' ~ 'Resistant')) %>% 
  ggplot(aes(y = categories, x = direction, fill = Category)) +
  geom_tile() +
  labs(
    x = 'Resistance',
    y = 'Terms'
  ) +
  scale_fill_manual(values = c('#2229F0', '#2292F0', '#22E4F0','#C7F0F0', 'white'),
                    name = 'FDR') +
  theme_classic()

ggsave(here('summary', 'aska_EcoCyc_enrichment_pval.pdf'), height =  5, width = 6.5)





# StringDB plots ----------------------------------------------------------

# load data
# bottom
bot_Function = read_delim("summary/StringDB_enrichment/bot.enrichment.Function.tsv", 
                                      "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
  rename(term = `term description`,
         fdr = `false discovery rate`) %>% 
  mutate(Direction = 'Bottom')

bot_KEGG = read_delim("summary/StringDB_enrichment/bot.enrichment.KEGG.tsv", 
                                  "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
  rename(term = `term description`,
         fdr = `false discovery rate`) %>% 
  mutate(Direction = 'Bottom')

bot_Net = read_delim("summary/StringDB_enrichment/bot.enrichment.NetworkNeighborAL.tsv", 
                      "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
  rename(term = `term description`,
         fdr = `false discovery rate`) %>% 
  mutate(Direction = 'Bottom')

bot_process = read_delim("summary/StringDB_enrichment/bot.enrichment.Process.tsv", 
                     "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
  rename(term = `term description`,
         fdr = `false discovery rate`) %>% 
  mutate(Direction = 'Bottom')


#top
top_Function = read_delim("summary/StringDB_enrichment/top.enrichment.Function.tsv", 
                          "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
  rename(term = `term description`,
         fdr = `false discovery rate`) %>% 
  mutate(Direction = 'Top')

top_KEGG = read_delim("summary/StringDB_enrichment/top.enrichment.KEGG.tsv", 
                      "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
  rename(term = `term description`,
         fdr = `false discovery rate`) %>% 
  mutate(Direction = 'Top')

top_Net = read_delim("summary/StringDB_enrichment/top.enrichment.NetworkNeighborAL.tsv", 
                     "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
  rename(term = `term description`,
         fdr = `false discovery rate`) %>% 
  mutate(Direction = 'Top')

top_process = read_delim("summary/StringDB_enrichment/top.enrichment.Process.tsv", 
                         "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
  rename(term = `term description`,
         fdr = `false discovery rate`) %>% 
  mutate(Direction = 'Top')










