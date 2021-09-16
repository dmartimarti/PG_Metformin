library(readr)
library(tidyverse)
library(here)
library(cowplot)

theme_set(theme_classic() +
            theme(axis.text.x = element_text(size = 13, color = 'black'),
                  axis.text.y = element_text(size = 13, color = 'black'),
                  axis.title.x = element_text(face = "bold", size = 13, color = 'black'),
                  axis.title.y = element_text(face = "bold", size = 13, color = 'black')))
# Read data ---------------------------------------------------------------


muts =  read_csv("D:/MRC_Postdoc/Pangenomic/mutants_analysis/results/mutations_summary.csv")

# prep the datatable
muts = muts %>% 
  mutate(Strain = as.factor(Strain), 
         TYPE = as.factor(TYPE),
         FTYPE = as.factor(FTYPE),
         STRAND = as.factor(STRAND)) %>% 
  separate(EFFECT, into = c('Variant', 'Variant_nt','Variant_aa'), sep = ' ') %>% 
  mutate(Variant = as.factor(Variant)) %>% 
  mutate(GENE = case_when(PRODUCT == 'hypothetical protein' ~ 'hypothetical protein',
                          TRUE ~ GENE)) %>% 
  replace_na(list( GENE = 'NO GENE'))

muts %>% write_csv(here('exploration','mutations_summary_wider.csv'))



#  effects on M and P variants --------------------------------------------


mp_muts = muts %>%
  filter(Strain %in% c('M_strains', 'P_strains'))


mp_muts %>% 
  # remove synonymous variants, not informative
  filter(Variant != 'synonymous_variant') %>%
  group_by(Strain) %>% 
  count(GENE) %>% 
  ggplot(aes(x = fct_reorder(GENE,n), y = n, fill = GENE)) +
  geom_histogram(stat="identity")  +
  labs(x = 'Gene variant',
       y = 'Number of total elements') +
  facet_wrap(~Strain,
             ncol = 1) +
  geom_text(aes(y = n+(2), x = GENE, label = round(n,0))) +
  theme_half_open(12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) 

ggsave(here('exploration','mutations_MP_strains.pdf'), height = 8, width = 10)


mp_muts %>% 
  # remove synonymous variants, not informative
  filter(Variant != 'synonymous_variant') %>%
  group_by(Strain, Colony, GENE) %>%
  summarise(N = n()) %>% 
  ggplot(aes(x = fct_reorder(GENE, N), y = N, fill = GENE)) +
  geom_histogram(stat="identity")  +
  labs(x = 'Gene variant',
       y = 'Number of total elements') +
  facet_wrap(~Strain*Colony)+
  geom_text(aes(y = N+(1), x = GENE, label = round(N,0))) + 
  theme_half_open(12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(here('exploration','mutations_MP_strains_big.pdf'), height = 18, width = 25)

mp_muts %>% 
  distinct(Colony, Strain) %>% 
  group_by(Strain) %>% 
  summary(N = n())

# calculate proportion of mutations per colony
# is a quick way to get a sense of how often they appear
mp_muts %>% 
  # remove synonymous variants, not informative
  filter(Variant != 'synonymous_variant') %>%
  group_by(Strain, GENE) %>%
  summarise(N = n()) %>% 
  group_by(Strain) %>% 
  # 19 P colonies, 12 M colonies
  mutate(prop = case_when(Strain == 'P_strains' ~ N / 19,
                          Strain == 'M_strains' ~ N / 12)) %>% 
  # remove rubish
  # filter(!(GENE %in% c('NO GENE','insB6_1', 'insA6_1'))) %>% 
  ggplot(aes(x = fct_reorder(GENE, prop), y = prop, fill = GENE)) +
  geom_histogram(stat="identity")  +
  labs(x = 'Gene variant',
       y = 'Proportion per colony') +
  geom_text(aes(y = prop+(0.1), x = GENE, label = round(prop,2))) + 
  geom_hline(yintercept = 1) +
  facet_wrap(~Strain,
             ncol = 1)+
  theme_half_open(12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(here('exploration','mut_props_MP_strains.pdf'), height = 8, width = 10)



### position of mutations within genes ####

# all mutations from rcdA gene
mp_muts %>% 
  separate(AA_POS, into = c('mut_pos', 'prot_length'), sep = '/', 
           convert = T, remove = F) %>% 
  filter(GENE == 'rcdA') %>% 
  group_by(GENE, Variant) %>% 
  count(mut_pos) %>% 
  left_join(expand_grid(mut_pos = 1:178)) %>% 
  ggplot(aes(x = mut_pos, y = n, fill = Variant, group = Variant)) +
  geom_histogram(stat = 'identity', color = 'black', size = 0.1) +
  xlim(0,178) +
  # facet_wrap(~GENE, scales = 'free_y') + 
  labs(x = 'Number of mutations',
       y = 'Position in the protein') +
  theme(legend.position="top")

ggsave(here('exploration','rcdA_mutations_pos.pdf'), height = 6, width = 9)



# missense mutations from rcdA genes
mp_muts %>% 
  separate(AA_POS, into = c('mut_pos', 'prot_length'), sep = '/', 
           convert = T, remove = F) %>% 
  filter(GENE == 'rcdA', Variant == 'missense_variant') %>% 
  group_by(GENE, Variant) %>% 
  count(mut_pos) %>% 
  left_join(expand_grid(mut_pos = 1:178)) %>% 
  ggplot(aes(x = mut_pos, y = n)) +
  geom_histogram(stat = 'identity', color = 'black', size = 0.4,
                 fill = 'dodgerblue3') +
  xlim(0,178) +
  #facet_wrap(~GENE, scales = 'free_y') +
  labs(x = 'Number of SNPs',
       y = 'Position in the protein') 

ggsave(here('exploration','rcdA_missense_mutations_pos.pdf'), height = 6, width = 9)





### combined plot with domains ####

doms = tibble('Domain'=c('DNA binding site', 'tetR domain'),
              'begin_pos'=c(30,13),
              'end_pos'=c(50,59))

muts_plot = mp_muts %>% 
  separate(AA_POS, into = c('mut_pos', 'prot_length'), sep = '/', 
           convert = T, remove = F) %>% 
  filter(GENE == 'rcdA', Variant == 'missense_variant') %>% 
  group_by(GENE, Variant) %>% 
  count(mut_pos) %>% 
  left_join(expand_grid(mut_pos = 1:178)) %>% 
  ggplot(aes(x = mut_pos, y = n)) +
  geom_histogram(stat = 'identity', color = 'black', size = 0.4,
                 fill = 'dodgerblue3') +
  xlim(0,178) +
  #facet_wrap(~GENE, scales = 'free_y') +
  labs(x = 'Number of SNPs',
       y = 'Position in the protein') +
  coord_cartesian(ylim = c(-0, 2), clip = "off")

doms_plot = ggplot(doms, aes(y = Domain)) +
  geom_segment(aes(x = begin_pos, xend = end_pos, yend = Domain,
                   color = Domain),
               size = 4) +
  guides(color='none') +
  annotate('text', x = 70, y = 2, label = 'tetR domain') +
  annotate('text', x = 70, y = 1, label = 'DNA binding') +
  xlim(0,178) +
  scale_color_manual(values = c('#F5D83C','#4B60F5')) +
  theme_void() +
  theme(
        axis.text.y=element_blank())

plot_grid(muts_plot, doms_plot, align = 'v', 
          rel_heights = c(10, 1),ncol = 1)

ggsave(here('exploration','rcdA_missense_mutations_pos_Domains.pdf'), height = 7, width = 9)


### 3D representation of PDB files


library(bio3d)

pdb = read.pdb(here('pdb_files', 'pdb7aco.ent'))

print(pdb)


modes <- nma(pdb)
plot(modes)

plot.bio3d(pdb$atom$b[pdb$calpha], 
           sse=pdb, typ="l", ylab="B-factor"
           )

# select only chain A
inds <- atom.select(pdb, chain=c("A"))

pdb2 <- trim.pdb(pdb, inds)


plot.bio3d(pdb2$atom$b[pdb2$calpha], 
           sse=pdb2, typ="l", ylab="B-factor"
)



rcdA_m1 = read.pdb(here('pdb_files', 'rcdA_M1_Gly51Glu.pdb'))
rcdA_wt = read.pdb(here('pdb_files', 'rcdA_WT.pdb'))

rcdA_aln = pdbaln(here('pdb_files', c('rcdA_M1_Gly51Glu.pdb',
                           'rcdA_WT.pdb')), fit = TRUE,
       web.args=list(email='dmartimarti@gmail.com'))



str_aln = struct.aln(rcdA_m1,rcdA_wt,
           web.args=list(email='dmartimarti@gmail.com'))


rcdA_m1.ind <- atom.select(rcdA_m1, chain="A", resno=15:55, elety="CA")
rcdA_wt.ind <- atom.select(rcdA_wt, chain="A", resno=15:55, elety="CA")

# perform superposition
xyz <- fit.xyz(fixed=rcdA_m1$xyz, mobile=rcdA_wt$xyz,
               fixed.inds=rcdA_m1.ind$xyz,
               mobile.inds=rcdA_wt.ind$xyz)

# write coordinates to file
write.pdb(rcdA_wt, xyz=xyz, file="rcdA_differences.pdb")


pdbs = read.pdb(here('pdb_files', c('rcdA_M1_Gly51Glu.pdb',
                                    'rcdA_WT.pdb')))

geostas(pdb)






#  effects on pangenome strains --------------------------------------------


pan_muts = muts %>%
  filter(!(Strain %in% c('M_strains', 'P_strains', 'Nissle')))



pan_muts %>% 
  # remove synonymous variants, not informative
  filter(Variant != 'synonymous_variant') %>%
  group_by(Strain) %>% 
  count(GENE) %>% 
  ggplot(aes(x = fct_reorder(GENE,n), y = n, fill = GENE)) +
  geom_histogram(stat="identity")  +
  labs(x = 'Gene variant',
       y = 'Number of total elements') +
  facet_wrap(~Strain,
             ncol = 3,
             scales = 'free_x') +
  geom_text(aes(y = n+(2), x = GENE, label = round(n,0))) + 
  theme_half_open(12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(here('exploration','mutations_Pangenome_strains.pdf'), height = 8, width = 10)



pan_muts %>% 
  # remove synonymous variants, not informative
  filter(Variant != 'synonymous_variant') %>%
  # group_by(Strain) %>% 
  count(GENE) %>% 
  ggplot(aes(x = fct_reorder(GENE,n), y = n, fill = GENE)) +
  geom_histogram(stat="identity")  +
  labs(x = 'Gene variant',
       y = 'Number of total elements') +
  # facet_wrap(~Strain,
  #            ncol = 3,
  #            scales = 'free_x') +
  geom_text(aes(y = n+(1), x = GENE, label = round(n,0))) + 
  theme_half_open(12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(here('exploration','mutations_Pangenome_total.pdf'), height = 8, width = 11)














