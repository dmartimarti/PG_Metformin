
# libraries ---------------------------------------------------------------


library(tidyverse)
library(readr)
library(broom)
library(gtools)
library(openxlsx)
library(here)
library(ComplexHeatmap)
library(rstatix)

#this sets up the plot theme, I hate the default grey background from ggplot
theme_set(theme_classic())

# read data

data = read_csv("Summary.csv") %>% 
  mutate(Type = factor(Type, 
                       levels = c('Treatment', 'Control')), # specify factor order
         Metabolite = as.factor(Metabolite),
         MetaboliteU = as.factor(MetaboliteU)) %>% 
  select(-File:-Data,-Comment,
         -`750nm_f_logAUC`:-`750nm_dt_Max`) %>% # remove shit 
  rename(AUC = `750nm_f_AUC`) # rename variable for easy use



# lm test -----------------------------------------------------------------

# a short test to demonstrate how linear model works
# filter the data to only one unique metabolite
test = data %>% filter(MetaboliteU == 'Potassium chromate|4')

# build the model
# dependent variable ~ independent variable
model = lm(AUC ~ Type, data = test)

# see how it worked
summary(model)

# the p-val for TypeControl is significant, so let's see in a boxplot

test %>% 
  ggplot(aes(x = Type, y = AUC, fill = Type)) + 
  geom_boxplot() +
  geom_point(position = position_jitterdodge()) # position jitterdodge will 
                                                # put every point in their 
                                                # groups, and not all together


# so it seems it works!



# multi-univariate statistics ---------------------------------------------

# let's work now with the entire dataset

# first, check that we have at least 2 samples per index
data %>% group_by(Index,Type) %>% count() %>% arrange(n)


# we are good to go, we have at least 2 samples per type and per metabolite

### drug conc diff ####

# this piece of code works in a very simple way:
# take the dataset, group by unique metabolite, nest the data into blocks, and 
# calculate the t-test per group
# unnest everything, tidy everything a bit, and just enjoy your new stats :)


stats = data %>% 
  group_by(MetaboliteU, Index) %>%                               # group by metabolites
  nest() %>%                                              # nest all data by the grouping variable
  mutate(lin_mod = map(data, lm, formula = 'AUC ~ Type'), # the map function will apply the function (lm) to the small datasets inside each row, see how the formula is input here
         lin_tidy = map(lin_mod, tidy)                    # this function (from the broom package) will tidy your results to look like a table
         )%>% 
  select(MetaboliteU, lin_tidy) %>%                       # select only the last variable you created, and the grouping variable
  unnest(lin_tidy) %>%                                    # unnest everything
  filter(term != '(Intercept)') %>%                       # remove this from the table, it's useless for us
  ungroup %>%                                             # ungroup the first group_by or the FDR won't work
  mutate(p.stars = stars.pval(p.value),                   # create a variable with stars depending on the p.value
         FDR = p.adjust(p.value, method = 'fdr'),         # calculate FDR
         FDR.stars = stars.pval(FDR))                     # new variable with stars but for the FDR values


# I did the comparison as Control vs Treatment, so get the inverse of 
# the estimate as a true Treament vs Control
stats = stats %>% 
  mutate(estimate = estimate * -1)

# look that everything is ok
stats


### total sum dif ####

# first we need to colapse the values into only one, by doing a sum of them
# (other methods can be tried later)

data.zip = data %>% 
  group_by(Metabolite, Type, Replicate) %>% 
  mutate(Total_sum = sum(AUC)) %>% 
  arrange(Replicate,MetaboliteU,Type) %>% 
  distinct(Total_sum, .keep_all = TRUE)


# do the stats
stats_sum = data.zip %>%
  group_by(Metabolite) %>%                               # group by metabolites
  nest() %>%                                              # nest all data by the grouping variable
  mutate(lin_mod = map(data, lm, formula = 'Total_sum ~ Type'), # the map function will apply the function (lm) to the small datasets inside each row, see how the formula is input here
         lin_tidy = map(lin_mod, tidy)                    # this function (from the broom package) will tidy your results to look like a table
  )%>% 
  select(Metabolite, lin_tidy) %>%                       # select only the last variable you created, and the grouping variable
  unnest(lin_tidy) %>%                                    # unnest everything
  filter(term != '(Intercept)') %>%                       # remove this from the table, it's useless for us
  ungroup %>%                                             # ungroup the first group_by or the FDR won't work
  mutate(p.stars = stars.pval(p.value),                   # create a variable with stars depending on the p.value
         FDR = p.adjust(p.value, method = 'fdr'),         # calculate FDR
         FDR.stars = stars.pval(FDR)) 

# I did the comparison as Control vs Treatment, so get the inverse of 
# the estimate as a true Treament vs Control
stats_sum = stats_sum %>% 
  mutate(estimate = estimate * -1)

# save everything into an excel file
list_of_tables = list(
  stats_R = stats,
  stats_sum = stats_sum
)

write.xlsx(list_of_tables, here('summary','Multi_univariate_stats.xlsx'))





# rstatix version ---------------------------------------------------------


stats = data %>% 
  group_by(MetaboliteU, Index) %>% 
  t_test(AUC ~ Type, detailed = TRUE) %>% 
  adjust_pvalue(method = "holm") %>%
  add_significance("p.adj") %>% 
  mutate(p.adj.stars = gtools::stars.pval(p.adj))



data.sum %>% 
  select(MetaboliteU, Type, Mean, SD) %>% 
  pivot_wider(values_from = c(Mean, SD), names_from = Type) %>% 
  ggplot(aes(x = Mean_Control, y = Mean_Treatment)) +
  geom_smooth(method = 'lm') + 
  geom_point(alpha = 0.5) +
  geom_errorbar(aes(ymax = Mean_Treatment + SD_Treatment, 
                    ymin = Mean_Treatment - SD_Treatment),
                alpha = 0.1) +
  geom_errorbarh(aes(xmax = Mean_Control + SD_Control, 
                     xmin = Mean_Control - SD_Control),
                 alpha = 0.1) +
  theme_classic()




# PCA ---------------------------------------------------------------------

### test with this data ####

data

library(tidymodels)

# create df for pca
data_pca = data %>% 
  unite(ID, Strain, Type, Replicate, sep = '_') %>% 
  unite(MetaboliteUPG, MetaboliteU, Index) %>% 
  select(ID, MetaboliteUPG, AUC) %>% 
  pivot_wider(names_from = MetaboliteUPG, values_from = AUC) 

# recipe and pca steps
rec <- recipe( ~ ., data = data_pca)
pca_trans <- rec %>%
  step_center(all_numeric()) %>%
  step_scale(all_numeric()) %>%
  step_pca(all_numeric(), num_comp = 3)
pca_estimates <- prep(pca_trans, training = data_pca)
pca_data <- bake(pca_estimates, data_pca)

# plot pca
pca_data %>% 
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_point() +
  geom_text(aes(label = ID),check_overlap = TRUE)

# % of var explained
sdev = pca_estimates$steps[[3]]$res$sdev
percent_variation = sdev^2 / sum(sdev^2)


### join full data ####

data_metf = read_csv("D:/MRC_Postdoc/Pangenomic/biolog/drug_drug_screen/summary/data_clean.csv")



# data wrangling to merge both datasets
data_metf_pca = data_metf %>% 
  filter(Plate != 'PM10') %>% 
  mutate(Strain = 'OP50',
         exp_number = 1) %>% 
  unite(MetaboliteUPG, MetaboliteU, Index) %>% 
  select(SampleID, MetaboliteUPG, AUC_raw) %>% 
  pivot_wider(names_from = MetaboliteUPG, values_from = AUC_raw) %>% 
  drop_na(`Cloxacillin|1_PM11C-B5`)

data_pca = data %>% 
  unite(ID, Strain, Type, Replicate, sep = '_') %>% 
  unite(MetaboliteUPG, MetaboliteU, Index) %>% 
  select(ID, MetaboliteUPG, AUC) %>% 
  pivot_wider(names_from = MetaboliteUPG, values_from = AUC) 

# datasets are not equal
dim(data_metf_pca)

dim(data_pca)

length(unique(data_metf$MetaboliteU))

length(unique(data$MetaboliteU))

#drop diff metabolites
drop_mets = sort(setdiff(data_metf$MetaboliteU, data$MetaboliteU))

data_metf_pca = data_metf %>% 
  filter(Plate != 'PM10') %>%
  filter(!(MetaboliteU %in% drop_mets)) %>% 
  mutate(Strain = 'OP50',
         exp_number = 1) %>% 
  unite(MetaboliteUPG, MetaboliteU, Index) %>% 
  select(SampleID, MetaboliteUPG, AUC_raw) %>% 
  pivot_wider(names_from = MetaboliteUPG, values_from = AUC_raw) %>% 
  drop_na(`Cloxacillin|1_PM11C-B5`)


dim(data_metf_pca)

dim(data_pca)

full_pca = rbind(data_metf_pca %>% select(ID = SampleID, everything()),
                 data_pca)
  

full_pca = full_pca %>% filter(ID != 'OP50_T_4_40') 

# recipe and pca steps
rec <- recipe( ~ ., data = full_pca)
pca_trans <- rec %>%
  step_center(all_numeric()) %>%
  step_scale(all_numeric()) %>%
  step_pca(all_numeric(), num_comp = 3)
pca_estimates <- prep(pca_trans, training = full_pca)
pca_data <- bake(pca_estimates, full_pca)




pca_data = pca_data %>% 
  mutate(Sample = case_when(
    str_detect(ID, "_60") ~ "Metf 60",
    str_detect(ID, "_40") ~ "Metf 40",
    str_detect(ID, "_20") ~ "Metf 20",
    str_detect(ID, "_0") ~ "Original Control",
    str_detect(ID, "OP50_Treatment") ~ "rcdA",
    str_detect(ID, "OP50_Control") ~ "New Control"
  ),
  Sample = as.factor(Sample))

# plot pca
pca_data %>% 
  ggplot(aes(x = PC1, y = PC2, color = Sample)) + 
  geom_point(size = 5)


ggsave(here('summary', 'pca_classes.pdf'), height = 8,width = 9)


# % of var explained
sdev = pca_estimates$steps[[3]]$res$sdev
percent_variation = sdev^2 / sum(sdev^2)

percent_variation[1:2]



# umap
library(embed)
umap_rec <- recipe(~., data = full_pca) %>%
  step_center(all_numeric()) %>%
  step_scale(all_numeric()) %>% 
  step_umap(all_numeric())

umap_prep <- prep(umap_rec)

umap_prep

umap_data = juice(umap_prep) %>% 
  mutate(Sample = case_when(
    str_detect(ID, "_60") ~ "Metf 60",
    str_detect(ID, "_40") ~ "Metf 40",
    str_detect(ID, "_20") ~ "Metf 20",
    str_detect(ID, "_0") ~ "Original Control",
    str_detect(ID, "OP50_Treatment") ~ "rcdA",
    str_detect(ID, "OP50_Control") ~ "New Control"
  ),
  Sample = as.factor(Sample))

umap_data %>%
  ggplot(aes(umap_1, umap_2, label = ID, color = Sample)) +
  # geom_bernie(bernie = 'sitting') +
  geom_point(size = 5) +
  # geom_text(check_overlap = TRUE) +
  labs(color = NULL)


ggsave(here('summary', 'umap_classes.pdf'), height = 8,width = 9)





# heatmap -----------------------------------------------------------------



# remove pm10 as it doesn't have MetaboliteU names for some of the compounds
data.sum = data %>% 
  group_by(Strain, Type, Index, MetaboliteU, Metabolite) %>%
  summarise(Mean = mean(AUC, na.rm = TRUE),
            SD = sd(AUC, na.rm = TRUE)) %>% 
  arrange(Type,MetaboliteU)


data.sum %>% 
  mutate(DrugConc = str_sub(MetaboliteU, -1))


### Heatmap ####

# lets try and create a heatmap with z-scores by metabolite and plate
library(ComplexHeatmap)

met.wide=data.sum %>%
  ungroup %>% 
  mutate(DrugConc = str_sub(MetaboliteU, -1),
         Plate = str_sub(Index, 1,4)) %>% 
  unite(MetID, Metabolite, Plate) %>%
  unite(ID, Strain, Type, DrugConc, sep = '_') %>% 
  select(ID, MetID, Mean) %>%
  pivot_wider(names_from = ID, values_from = Mean) %>% 
  select(MetID, OP50_Control_1:OP50_Control_4,OP50_Treatment_1:OP50_Treatment_4)

# some stupid data formatting, change met.wide from tibble to data.frame
# the objective here is to have a matrix with row names as metabolites, 
# and col names as samples
# IMPORTANT: matrix is NOT a data frame, and data frame is NOT a tibble (format for tidyverse functions)
met.wide = data.frame(met.wide)
rownames(met.wide) = met.wide[,1] # set row names as metabolites
met.wide[,1] = NULL # remove metabolite column, 

# replace Inf values for NAs, just in case 
met.wide[sapply(met.wide, is.infinite)] <- NA

# scale values! this is how you get the z-score
# here I'm using the transpose of the matrix because if not, we would get the z-score by samples (columns) instead
met.wide.scale = scale(t(met.wide))

# annotation layer
ha = columnAnnotation(strain = c(rep('OP50',4),rep('rcdA',4)),
                      drug_conc = c(1,2,3,4,1,2,3,4),
                      col = list(strain = c("rcdA" = "red", "OP50" = "green")),
                      border = TRUE)

# plot with Complexheatmaps library, it can be done with ggplot as well
Heatmap(t(met.wide.scale), # again, transpose the matrix to represent data in a way we can understand
        name = "Z-score", # name of the legend
        column_title = "Samples",
        row_title = "Metabolites",
        row_names_gp = gpar(fontsize = 5),
        cluster_columns = FALSE,
        top_annotation = ha) # reduce font size in rows






# save plot
dev.copy2pdf(device = cairo_pdf,
             file = 'Heatmap.pdf',
             width = 10, height = 15, useDingbats = FALSE)




# boxplots ----------------------------------------------------------------

drug = 'Spectinomycin'

# tests
data %>%
  mutate(DrugConc = str_sub(MetaboliteU, -1),
         Plate = str_sub(Index, 1,4),
         Replicate = as.factor(Replicate)) %>%
  filter(Metabolite == drug) %>% 
  mutate(Type = factor(Type, levels = c('Control', 'Treatment'))) %>% 
  ggplot(aes(x = DrugConc, y = AUC, fill = Type)) + 
  geom_boxplot() +
  geom_point(position = position_jitterdodge()) 

data.sum %>%
  filter(Metabolite == drug) %>% 
  mutate(Type = factor(Type, levels = c('Control', 'Treatment'))) %>% 
  ggplot(aes(x = MetaboliteU, y = Mean, color = Type, fill = Type, group = Type)) + 
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.1) +
  geom_line() + geom_point()

# BIG PLOTS

p = data %>%
  mutate(DrugConc = str_sub(MetaboliteU, -1),
         Plate = str_sub(Index, 1,4),
         Replicate = as.factor(Replicate)) %>% 
  unite(MetIndex, Metabolite, Plate, remove=F) %>% 
  mutate(Type = factor(Type, levels = c('Control', 'Treatment'))) %>% 
  ggplot(aes(x = DrugConc, y = AUC, fill = Type)) + 
  geom_boxplot() +
  geom_point(aes(color = Replicate), position = position_jitterdodge()) +
  facet_wrap(~MetIndex, scales = 'free') 


ggsave(here('summary','boxplots_COMPLETE.pdf'), height = 25, width = 35)



data.sum %>%
  mutate(DrugConc = str_sub(MetaboliteU, -1),
         Plate = str_sub(Index, 1,4)) %>% 
  unite(MetIndex, Metabolite, Plate, remove=F) %>% 
  mutate(Type = factor(Type, levels = c('Control', 'Treatment'))) %>% 
  ggplot(aes(x = DrugConc, y = Mean, color = Type, fill = Type, group = Type)) + 
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.1) +
  geom_line() + 
  geom_point() + 
  facet_wrap(~MetIndex, scales = 'free')

ggsave(here('summary','LinePlot_COMPLETE.pdf'), height = 25, width = 35)





# enrichment --------------------------------------------------------------



library(readxl)
biolog = read_excel("D:/MRC_Postdoc/Pangenomic/biolog/drug_drug_screen/Biolog_metabolites_UPDATED_drugClass.xlsx", 
                                                   sheet = "Biolog_drugs")


# make a table with drug classes for enrichment data
drug_target = biolog %>%
  mutate(DrugConc = str_sub(MetaboliteU, -1)) %>%
  filter(DrugConc == 1) %>%
  select(Plate, Well, Index, Metabolite, EcoCycID, KEGG_ID, Target) %>%
  separate_rows(Target, sep = ', ') %>%
  drop_na(Target) %>%
  mutate(Target = as.factor(Target)) %>%
  unite(Drug.combination, Metabolite, Plate, remove = F)

data %>% 
  mutate(DrugConc = str_sub(MetaboliteU, -1),
         Plate = str_sub(Index, 1,4))




enrich = function(syn.items=syn, ant.items=ant, db=biolog, feature='Target'){
  # initialise variables
  drug_target = biolog %>%
    mutate(DrugConc = str_sub(MetaboliteU, -1)) %>%
    filter(DrugConc == 1) %>%
    select(Plate, Well, Index, Metabolite, EcoCycID, KEGG_ID, feature) %>%
    separate_rows(feature, sep = ', ') %>%
    drop_na(feature) %>%
    unite(Drug.combination, Metabolite, Plate, remove = F)
  
  
  classes = drug_target %>% select(feature) %>% t %>% as.vector %>% unique
  N = length(unique(drug_target$Metabolite))
  
  # hypergeometric test
  # synergy
  syn.enrich = c()
  for (class in classes){
    class.met = drug_target %>% filter(!!as.symbol(feature) == class) %>% select(Drug.combination) %>% t %>% as.vector
    m = length(class.met)
    n = N - m
    k = length(syn.items)
    x = length(class.met[class.met %in% syn.items])
    fit = phyper(q = x-1, m = m, n = n, k = k, lower.tail = FALSE)
    syn.enrich = c(syn.enrich, fit)
  }
  
  # antagonistic
  ant.enrich = c()
  for (class in classes){
    class.met = drug_target %>% filter(!!as.symbol(feature) == class) %>% select(Drug.combination) %>% t %>% as.vector
    m = length(class.met)
    n = N - m
    k = length(ant.items)
    x = length(class.met[class.met %in% ant.items])
    fit = phyper(q = x-1, m = m, n = n, k = k, lower.tail = FALSE)
    ant.enrich = c(ant.enrich, fit)
  }
  
  df = data.frame(classes,syn.enrich, ant.enrich)
  colnames(df) = c('Class', 'Synergy', 'Antagonism')
  df = df %>% mutate(Syn.stars = stars.pval(Synergy),
                     Ant.stars = stars.pval(Antagonism)) %>%
    select(Class, Synergy, Syn.stars, Antagonism, Ant.stars)
  return(tibble(df))
}


# divide data in synergistic/antagonistic
# IMPORTANT, READ!!
# as my script takes metabolite names + plate name, I need to create those

enrich_data = stats_sum %>% 
  left_join(biolog %>% 
              select(-MetaboliteU) %>% 
              distinct(Metabolite,.keep_all=TRUE)) %>% 
  unite(Metab, Metabolite, Plate) 



ant = enrich_data %>% 
  filter(estimate > 0, p.value < 0.05) %>% select(Metab) %>% t %>% 
  as.vector

syn = enrich_data %>% 
  filter(estimate < 0, p.value < 0.05) %>% select(Metab) %>% t %>% 
  as.vector

N = length(unique(data$Metabolite))


# calculate results
mol.res = enrich(syn, ant, biolog, 'Molecule')
target.res = enrich(syn, ant, biolog, 'Target') 
process.res = enrich(syn, ant, biolog, 'Process')







# plot
sig = 0.05
df = target.res %>%
  mutate(Comparison = 'Target') %>%
  bind_rows(mol.res %>% mutate(Comparison = 'Molecule'), 
            process.res %>% mutate(Comparison = 'Process')) %>%
  filter(Synergy < sig | Antagonism < sig) %>%
  select(Class, Synergy, Antagonism, Comparison)  %>%
  pivot_longer(cols = c('Synergy', 'Antagonism'), names_to = 'Direction', values_to = 'p.value')


# enrichment procedure
enrbrks = c(0, -log(0.1, 10), -log(0.05, 10), 2, 3, 4, 100)
enrlbls = c('N.S.', '<0.1', '<0.05','<0.01','<0.001','<0.0001')
enrcols = colorRampPalette(c("gray90", "steelblue1", "blue4"))(n = 7)


p.theme = theme(axis.ticks = element_blank(), panel.border = element_blank(), 
                panel.background = element_blank(), panel.grid.minor = element_blank(), 
                panel.grid.major = element_blank(), axis.line = element_line(colour = NA), 
                axis.line.x = element_line(colour = NA), axis.line.y = element_line(colour = NA), 
                strip.text = element_text(colour = "black", face = "bold", 
                                          size = 7), 
                axis.text.x = element_text(face = "bold", 
                                           colour = "black", size = 10, angle = 45, hjust = 1))


# plot enrichment p-values
df %>% 
  mutate(p.value = p.value + 0.00000001,
         logFDR = ifelse(-log10(p.value) < 0, 0, -log10(p.value)),
         logFDRbin = cut(logFDR, breaks = enrbrks, labels = enrlbls, right = FALSE),
         Class = factor(Class),
         Direction = factor(Direction)) %>%
  ggplot(aes(x = Direction, y = Class)) +
  geom_tile(aes(fill = logFDRbin)) +
  scale_fill_manual(values = enrcols)  + 
  facet_wrap(~Comparison, scales = 'free_y') +
  p.theme


dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'Drug_enrichment_relaxed.pdf'),
             width = 8, height = 7, useDingbats = FALSE)






# Growth curves (sig metab) -----------------------------------------------








# Get timeseries data
# It will take a while
data.t = read_csv('Timeseries.csv',quote = "\"") %>%
  filter(Data == '750nm_f') %>% 
  gather(Time_s, OD, matches('\\d')) %>%
  filter(!is.na(OD)) %>% 
  mutate(Type = ifelse(Type == 'Control','C','T'),
         Type = factor(Type,
                       levels = c('C','T'),
                       labels = c('Control','Treatment')),
         Time_s = as.numeric(Time_s),
         Time_h = Time_s/3600,
         Row = str_match_all(Well,'[:digit:]{1,}'),
         Col = str_match_all(Well,'[:alpha:]{1,}'), 
         Row = factor(Row, levels = 1:12), 
         Col = factor(Col, levels = LETTERS[1:8])) %>% 
  select(-c(File,Reader, Comment))

# Growth curves for summary

tsum = data.t %>%
  group_by(Strain, Type, Index, Plate, Well, Metabolite, MetaboliteU, Time_h) %>%
  summarise(Mean = mean(OD), SD = sd(OD),SE = SD/sqrt(length(OD))) %>%
  ungroup 


# test to plot a single metabolite
tsum %>%
  filter(Metabolite == 'Potassium chromate') %>% 
  ggplot( aes(x = Time_h, y = Mean, fill = Type, color = Type)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.2) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 24, by = 12)) +
  ylab("OD") +
  xlab("Time, h") +
  labs(fill = "Type") +
  facet_grid((~MetaboliteU*Metabolite)) 


# let's select only the significant metabolites
sig_metabs = stats_sum %>%
  filter(p.value < 0.05) %>% 
  select(Metabolite) %>% t %>% as.character


# This is a very busy and useless plot

# # plot all sig metabs growth curves
# tsum %>%
#   filter(Metabolite %in% sig_metabs) %>% 
#   ggplot( aes(x = Time_h, y = Mean, fill = Type, color = Type)) +
#   geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.2) +
#   geom_line() +
#   scale_x_continuous(breaks = seq(0, 24, by = 12)) +
#   ylab("OD") +
#   viridis::scale_fill_viridis(discrete=TRUE) +
#   viridis::scale_color_viridis(discrete=TRUE) +
#   xlab("Time, h") +
#   labs(fill = "Type") +
#   facet_wrap((~MetaboliteU), ncol = 4) 
# 
# 
# dev.copy2pdf(device = cairo_pdf,
#              file = here('summary', 'Growth_curves_sig_metab.pdf'),
#              width = 10, height = 60, useDingbats = FALSE)


# THIS LOOP PLOTS EACH METABOLITE IN A SEPARATE PLOT

for ( met in sig_metabs){

  # plot all sig metabs growth curves
  tsum %>%
    filter(Metabolite %in% met) %>% 
    ggplot( aes(x = Time_h, y = Mean, fill = Type, color = Type)) +
    geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.2) +
    geom_line() +
    scale_x_continuous(breaks = seq(0, 24, by = 12)) +
    ylab("OD") +
    viridis::scale_fill_viridis(discrete=TRUE) +
    viridis::scale_color_viridis(discrete=TRUE) +
    xlab("Time, h") +
    labs(fill = "Type") +
    facet_wrap((~MetaboliteU), ncol = 4) 
  
  
  ggsave(here('summary/growth_curves', paste0(met,'_growth_curves.pdf')),
               width = 10, height = 8)
}

### MetaboliteU filter

# let's select only the significant metabolites
sig_metabs = stats %>%
  filter(p.value < 0.05) %>% 
  select(MetaboliteU) %>% 
  mutate(MetaboliteU = str_sub(MetaboliteU,1,-3)) %>% 
  distinct(MetaboliteU) %>% 
  t %>% as.character


# THIS LOOP PLOTS EACH METABOLITE IN A SEPARATE PLOT

for ( met in sig_metabs){
  
  # plot all sig metabs growth curves
  tsum %>%
    filter(Metabolite %in% met) %>% 
    ggplot( aes(x = Time_h, y = Mean, fill = Type, color = Type)) +
    geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.2) +
    geom_line() +
    scale_x_continuous(breaks = seq(0, 24, by = 12)) +
    ylab("OD") +
    viridis::scale_fill_viridis(discrete=TRUE) +
    viridis::scale_color_viridis(discrete=TRUE) +
    xlab("Time, h") +
    labs(fill = "Type") +
    facet_wrap((~MetaboliteU), ncol = 4) 
  
  
  ggsave(here('summary/growth_curves_metabU', paste0(met,'_growth_curves.pdf')),
         width = 10, height = 8)
}







