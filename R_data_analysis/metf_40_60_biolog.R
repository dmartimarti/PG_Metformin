# libraries ####
library(tidyverse) # master library to deal with data frames
library(readxl) # read xlsx or xls files
library(ggrepel) # ggplot add-on, to plot names that don't collapse in same position
library(FactoMineR) # for PCA
library(factoextra) # for PCA
library(here) # usefull to save plots in folders given a root
library(viridis) # color palette package
library(ComplexHeatmap) # yeah, complex heatmaps
library(broom)
library(PFun)
# set wd
# I unpacked Sim’s biolog file here, you can do it wherever you want in your computers
setwd("D:/MRC_Postdoc/Pangenomic/biolog/Biolog_metf_40_60")

### Load data from sources ####

# load every file from the list of folders
dir = getwd()

plates = c('pm10', 'pm11', 'pm12', 'pm13', 'pm14','pm15',
           'pm16', 'pm17', 'pm18', 'pm19', 'pm20')

paths = file.path(dir, 'biologdata', c('40mM', '60mM'))

paths40 = paste(paths[1], plates, 'Output/Summary.csv', sep = '/')
paths60 = paste(paths[2], plates, 'Output/Summary.csv', sep = '/')

###
# get data from 40mM files
###


data40 = data.frame()
for (path in paths40){
  temp = read_csv(path, quote = "\"")
  data40 = rbind(data40, temp)
}

# put the data in shape
data40 = data40 %>% rename(AUC_raw = `750nm_f_AUC`, 
                       Metformin_mM = Comment) %>% # `750nm_f_logAUC` data column is what we need for logAUC values
  mutate(Strain = as.character(Strain),
         Type = ifelse(Type == 'Control','C','T'),
         Sample = paste(Strain, Type, sep = '_'),
         Type = factor(Type,
                       levels = c('C', 'T'),
                       labels = c('Control', 'Treatment')),
         Row = str_match_all(Well,'[:digit:]{1,}'), #Get plate row names from well name. 
         Col = str_match_all(Well,'[:alpha:]{1,}'), #Get plate column names from well name
         Row = factor(Row, levels = 1:12), #Make them categorical variable with set order them
         Col = factor(Col, levels = LETTERS[1:8]),
         Metformin_mM = ifelse(Metformin_mM == 'zero', 0, ifelse(Metformin_mM == '60mM', 60, 40)),
         SampleID = paste(Sample, Replicate, Metformin_mM, sep = '_'),
         Sample = as.factor(Sample),
         Strain = as.factor(Strain)) %>% #Change Type column coding
  select(SampleID, Metformin_mM, Sample, Strain, Type, Replicate, Replicate, Index, Plate, 
         Well, Row, Col, Metabolite, MetaboliteU, EcoCycID, KEGG_ID, Group, Class, AUC_raw)

###
# get data from 60mM files
###


data60 = data.frame()
for (path in paths60){
  temp = read_csv(path, quote = "\"")
  data60 = rbind(data60, temp)
}


# put the data in shape
data60 = data60 %>% rename(AUC_raw = `750nm_f_AUC`, 
                         Metformin_mM = Comment) %>% # `750nm_f_logAUC` data column is what we need for logAUC values
  mutate(Strain = as.character(Strain),
         Type = ifelse(Type == 'Control','C','T'),
         Sample = paste(Strain, Type, sep = '_'),
         Type = factor(Type,
                       levels = c('C', 'T'),
                       labels = c('Control', 'Treatment')),
         Row = str_match_all(Well,'[:digit:]{1,}'), #Get plate row names from well name. 
         Col = str_match_all(Well,'[:alpha:]{1,}'), #Get plate column names from well name
         Row = factor(Row, levels = 1:12), #Make them categorical variable with set order them
         Col = factor(Col, levels = LETTERS[1:8]),
         Metformin_mM = ifelse(Metformin_mM == 'zero', 0, ifelse(Metformin_mM == '60mM', 60, 40)),
         SampleID = paste(Sample, Replicate, Metformin_mM, sep = '_'),
         Sample = as.factor(Sample),
         Strain = as.factor(Strain)) %>% #Change Type column coding
  filter(Metformin_mM != 0) %>%
  select(SampleID, Metformin_mM, Sample, Strain, Type, Replicate, Replicate, Index, Plate, 
         Well, Row, Col, Metabolite, MetaboliteU, EcoCycID, KEGG_ID, Group, Class, AUC_raw)


# join both datasets
data = rbind(data40, data60)

# set a column with drug concentrations
data = data %>%
  mutate(DrugConc = str_sub(MetaboliteU, -1)) %>%
  unite(Type_met, Type, Metformin_mM, remove = FALSE)


### Summary stats ####


data.sum = data %>% 
  group_by(Metformin_mM, Strain, Sample, Type, Index, Plate, MetaboliteU, Metabolite) %>%
  summarise(Mean = mean(AUC_raw, na.rm = TRUE),
            SD = sd(AUC_raw, na.rm = TRUE))

cosa = data.sum %>% ungroup %>%
  select(Sample, Metformin_mM, MetaboliteU, Mean) %>%
  unite(ID, Sample, Metformin_mM)

mets = cosa %>% filter(ID == 'OP50_C_0') %>% select(MetaboliteU) %>% t %>% as.character
dup.mets = mets[duplicated(mets)]

data.sum.wide = data.sum %>% ungroup %>%
  filter(!MetaboliteU %in% dup.mets) %>%
  select(Sample, Metformin_mM, MetaboliteU, Mean) %>%
  unite(ID, Sample, Metformin_mM) %>%
  pivot_wider(names_from = ID, values_from = Mean)

mets.fix = data.sum.wide %>% select(MetaboliteU) %>% t %>% as.character

rescue40 = data.sum.wide[,3] - data.sum.wide[,2]
met40 = mets.fix[rescue40 >= 1]

rescue60 = data.sum.wide[,4] - data.sum.wide[,2]
met60 = mets.fix[rescue60 >= 1]

data %>%
  filter(MetaboliteU %in% met40) %>%
  filter(Metformin_mM %in% c(0,40)) %>%
  mutate(Metformin_mM = as.factor(Metformin_mM)) %>%
  ggplot(aes(x = MetaboliteU, y = AUC_raw, fill = Metformin_mM)) +
  geom_boxplot() +
  # geom_point() +
  theme_classic() +
  geom_vline(xintercept = c(1.5:length(met40)), color = "gray40", size = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.copy2pdf(device = cairo_pdf,
             file ='test.pdf',
             height = 8, width = 15, useDingbats = FALSE)

sub = data %>% 
  filter(MetaboliteU %in% met40[1]) %>%
  mutate(Metformin_mM = as.factor(Metformin_mM))
model = aov(AUC_raw ~ Metformin_mM, data = sub)
TukeyHSD(model)

tidy(TukeyHSD(model)) %>% filter(comparison == '40-0') %>% mutate(MetaboliteU = 'cosa')

cosa = data.frame()
for (met in met40){
  sub = data %>% filter(MetaboliteU == met) %>%
    mutate(Metformin_mM = as.factor(Metformin_mM))
  model = aov(AUC_raw ~ Metformin_mM, data = sub)
  temp = tidy(TukeyHSD(model)) %>% filter(comparison == '40-0') %>% mutate(MetaboliteU = met)
  cosa = rbind(cosa, temp)
}

cosa = cosa %>% mutate(adj.p.value = p.adjust(adj.p.value, method = 'fdr'))






# lets try and create a heatmap with z-scores by metabolite and plate
# this will imply some data shaping and manipulation

# first step: make data wide. Select uracil == 40 and remove Neg control
# select only the variables you want to keep, and make it wide
data.wide = data %>%
  filter(Plate != 'PM10') %>%
  filter(Metformin_mM %in% c(40)) %>%
  select(SampleID, MetaboliteU, AUC_raw) %>%
  pivot_wider(names_from = c(SampleID), values_from = AUC_raw)

# some stupid data formatting, change met.wide from tibble to data.frame
# the objective here is to have a matrix with row names as metabolites, 
# and col names as samples
# IMPORTANT: matrix is NOT a data frame, and data frame is NOT a tibble (format for tidyverse functions)
met.wide = data.frame(met.wide)
rownames(met.wide) = met.wide[,1] # set row names as metabolites
met.wide[,1] = NULL # remove metabolite column, 

# replace Inf values for NAs, need to check what happend here
met.wide[sapply(met.wide, is.infinite)] <- NA

# scale values! this is how you get the z-score
# here I'm using the transpose of the matrix because if not, we would get the z-score by samples (columns) instead
met.wide.scale = scale(t(met.wide))

# plot with Complexheatmaps library, it can be done with ggplot as well
Heatmap(t(met.wide.scale), # again, transpose the matrix to represent data in a way we can understand
        name = "Z-score", # name of the legend
        column_title = "Samples",
        row_title = "Metabolites",
        row_names_gp = gpar(fontsize = 5)) # reduce font size in rows

# save plot
dev.copy2pdf(device = cairo_pdf,
             file = 'Heatmap.pdf',
             width = 10, height = 15, useDingbats = FALSE)









# # # # # # # # # # # # # # # #  
### Multi-univariate stats ####
# # # # # # # # # # # # # # # # 


# let's try to calculate interactin terms with drugs and with emmeans package

library(emmeans)

sub = data %>% filter(Metabolite == 'D−Cycloserine') %>%
  dplyr::select(Metformin_mM, MetaboliteU, Metabolite, AUC_raw)

sub

sub.lm <- lm(AUC_raw ~ MetaboliteU * Metformin_mM, data = sub)
anova(sub.lm)

emmip(sub.lm, Metformin_mM~MetaboliteU)

emmeans(sub.lm, pairwise ~ Metformin_mM * MetaboliteU)


# check metabolite behaviour
met = 'D−Cycloserine'
stats %>% filter(Metabolite == met)
data %>% filter(Metabolite == met) %>%
  ggplot(aes(x = Type_met, y = AUC_raw, fill = Metformin_mM)) +
  geom_boxplot() +
  facet_wrap(~MetaboliteU) +
  theme_classic()










####
# Linear Mod with Pov's function

contrasts <- read.contrasts('!Contrasts.xlsx') # This table describes comparisons made in data

# loads the metadata table
contrasts.table <- contrasts$Contrasts.table
# loads the contrasts matrix (1s and 0s)
contr.matrix <- contrasts$Contrasts.matrix

# shortens the metadata data frame
contr.desc <- contrasts.table %>%
  select(Description:Strain)



contr.matrix

# Initialise adjustments
group <- c('All',"PM11C", "PM12B", "PM13B", "PM14A", "PM15B", "PM16A", "PM17A", "PM18C", "PM19" , "PM20B")
strain <- c('OP50')


adjustments <- expand.grid('Group' = group, 'Strain' = strain) %>%
  mutate(a = -1, b = 0,
         Contrast = paste0(Strain, '_Metf_', Group, 'n')) %>%
  select(Contrast, Group:b) 



# transform the contr.matrix into a real matrix (not data frame)
contr.matrix <- contr.matrix %>%
  as.matrix()


# this piece of code performs a lot of different things, be careful

allresults <- data %>%
  unite(Met_ID, Metabolite, Plate, remove = FALSE) %>%
  #filter(Strain == 'BW') %>%
  filter(!Metabolite %in% c('Negative Control','Positive Control') ) %>%
  # filter(! (Plate == 'PM3B' & SampleID == 'BW_C_1') ) %>%
  # filter(! (Plate == 'PM4A' & SampleID == 'TM_C_3') ) %>%
  group_by(Plate, Well, Index, Met_ID, Metabolite, MetaboliteU, EcoCycID, KEGG_ID) %>%
  # Another creation of mine. Nothing fancy, just a wrapper for the use of contrasts in LM. 
  # But integrates nicely with tidyverse workflow. Code in PFun
  do(hypothesise(.,"AUC_raw~0+Type_met", contr.matrix)) %>% 
  group_by(Plate, Well, Met_ID, Metabolite, MetaboliteU, EcoCycID, KEGG_ID) %>%
  getresults(contr.desc) # Calculate some additional parameters. Code in PFun

#get results in different shapes
results <- allresults$results
results.cast <- allresults$cast
results.castfull <- allresults$castfull

results$Contrast %>% unique %>% as.character


results %>% filter(Contrast_type == 'Strain')

# Save statistical analysis results
write.csv(results,here('summary','Metf_results.csv'),row.names = FALSE)
write.csv(results.cast,here('summary','Metf_results_sidebyside.csv', sep = ''), row.names = FALSE)
write.csv(results.castfull,here('summary','Metf_results_sidebyside_full.csv',sep=''), row.names = FALSE)



### line Plot ####

data.sum %>%
  ungroup %>%
  mutate(Metformin_mM = as.factor(Metformin_mM),
         DrugConc = str_sub(MetaboliteU, -1)) %>%
  filter(Plate == 'PM11C') %>%
  ggplot(aes(x = DrugConc, y = Mean, colour = Metformin_mM, fill = Metformin_mM, group = Metformin_mM)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), color = NA, alpha = 0.2) +
  geom_line(size = 1) +
  theme_classic() +
  facet_wrap(~Metabolite) +
  theme(strip.text = element_text(size = 5),
        plot.title = element_text(hjust = 0.5, size = 5,
                                  color = 'black')) 

# # # # # # # # # # # # # # 
### SynergyFinder test ####
# # # # # # # # # # # # # #

# Erythromycin, D-Cycloserine: antagonistic
# Atropine, Orphenadrine; synergistic
# Amikacin: neutral
drugs = c('Amikacin', 'D-Cycloserine', 'Erythromycin', 'Atropine', 'Orphenadrine')
other_ctr = c('X-b-D-Glucuronide', 'X-Caprylate')
pH = c('pH 5.5', 'pH 6', 'pH 6.5', 'pH 7')

data %>% filter(Plate == 'PM10', Metabolite %in% other_ctr) %>%
  mutate(Metformin_mM = as.factor(Metformin_mM),
         Metabolite = as.factor(Metabolite)) %>%
  ggplot(aes(x = Metabolite, y = AUC_raw, fill = Metformin_mM)) +
  geom_boxplot() 


# let's use pH7 as control for drug == 0, metf == (0,40,60)

# substracting drugs and renaming everything
test = data %>% filter(Metabolite %in% drugs) %>%
  mutate(DrugConc = str_sub(MetaboliteU, -1),
         Drug1 = 'Metformin',
         ConcUnit = 'A.U.',
         PairIndex = 1) %>%
  select(PairIndex, Metformin_mM, Metabolite, DrugConc, AUC_raw, Drug1, ConcUnit) %>%
  rename(Conc1 = Metformin_mM,
         Conc2 = DrugConc, 
         Drug2 = Metabolite,
         Response = AUC_raw) %>%
  arrange(Drug2, Conc1, Conc2) %>%
  mutate(Conc2 = as.numeric(Conc2)) %>%
  select(PairIndex, Drug1, Drug2, Conc1, Conc2, Response, ConcUnit)

# similar treatment to get ph 7 as controls
drug_ctr = data %>% filter(Metabolite %in% c('X-b-D-Glucuronide')) %>%
  mutate(DrugConc = str_sub(MetaboliteU, -1),
         DrugConc = 0,
         Drug1 = 'Metformin',
         ConcUnit = 'A.U.',
         PairIndex = 1) %>%
  select(PairIndex, Metformin_mM, Metabolite, DrugConc, AUC_raw, Drug1, ConcUnit) %>%
  rename(Conc1 = Metformin_mM,
         Conc2 = DrugConc,
         Drug2 = Metabolite,
         Response = AUC_raw) %>%
  arrange(Drug2, Conc1, Conc2) %>%
  select(PairIndex, Drug1, Conc1, Conc2, Response, ConcUnit) 

cosa = expand.grid(Drug1 = 'Metformin', Drug2 = drugs)

complete.test = cosa %>% left_join(drug_ctr) %>% rbind(test) %>% tibble

# take the mean of pH 6 to generate the % of inhibition
mean_ctrl = data %>% filter(Metabolite == 'pH 6', Metformin_mM == 0) %>%
  summarise(Mean = mean(AUC_raw))

# calculate % of inhibition, and remake table to save into csv
complete.test = complete.test %>% 
  mutate(inhibition = 100 - (Response/14.4 * 100)) %>%
  select(-Response) %>%
  rename(Response = inhibition) %>%
  mutate(PairIndex = 1) %>%
  arrange(Drug2)


write.csv(complete.test, 'SynergyFinder_test.csv', row.names = FALSE)

complete.test2 = complete.test %>% mutate(Conc1 = ifelse(Conc1 == 40, 4, ifelse(Conc1 == 60, 6, 0)))
write.csv(complete.test2, 'SynergyFinder_test_conc.csv', row.names = FALSE)


# second 

data %>%
  filter(MetaboliteU %in% met40) %>%
  filter(Metformin_mM %in% c(0,40)) %>%
  mutate(Metformin_mM = as.factor(Metformin_mM)) %>%
  ggplot(aes(x = MetaboliteU, y = AUC_raw, fill = Metformin_mM)) +
  geom_boxplot() +
  # geom_point() +
  theme_classic() +
  geom_vline(xintercept = c(1.5:length(met40)), color = "gray40", size = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


### exploration to select the best control ####
# - plot all drugs at concentration == 1
# - select the ones that have the higher AUC
# - see how they behave with metf 40 and 60
# - decide which one to Use
# - establish a pipeline to automate this thing



data %>% filter(DrugConc == 1, Metformin_mM == 0) %>%
  group_by(Metabolite) %>%
  summarise(Mean = mean(AUC_raw),
            SD = sd(AUC_raw),
            SEM = SD/sqrt(n())) %>%
  arrange(desc(Mean, SD))
















