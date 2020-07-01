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

### custom functions ####

# I extracted this one from the synergyfinder package
ExtendedScores <-  function (scores.mat, len) {
  # len: how many values need to be predicted between two adjacent elements
  #      of scores.dose
  options(scipen = 999)
  nr <- nrow(scores.mat)
  nc <- ncol(scores.mat)
  
  # missing value imputation
  while (sum(is.na(scores.mat))) {
    scores.mat <- ImputeNA(scores.mat)
  }
  ext.row.len <- (nr - 1) * (len + 2) - (nr - 2)
  ext.col.len <- (nc - 1) * (len + 2) - (nc - 2)
  
  extended.row.idx <- seq(1, nr, length = ext.row.len)
  extended.col.idx <- seq(1, nc, length = ext.col.len)
  
  krig.coord <- cbind(rep(extended.row.idx, each = ext.col.len),
                      rep(extended.col.idx, times = ext.row.len))
  extended.scores <- kriging(data = c(scores.mat),
                             data.coord = cbind(rep(seq_len(nr), nc),
                                                rep(seq_len(nc), each=nr)),
                             krig.coord = krig.coord,
                             cov.mod = "whitmat", grid = FALSE,
                             sill = 1, range = 10,
                             smooth = 0.8)$krig.est
  extended.scores <- matrix(extended.scores, nrow = ext.row.len,
                            ncol = ext.col.len, byrow = TRUE)
  
  # extended.scores = data.frame(extended.scores)
  extended.scores <- round(extended.scores, 3)
  
  row.dose <- as.numeric(rownames(scores.mat))
  col.dose <- as.numeric(colnames(scores.mat))
  
  extend.row.dose <- mapply(function(x, y){seq(from = x, to = y,
                                               length.out = len + 2)},
                            row.dose[-nr], row.dose[-1])
  extend.row.dose <- unique(round(c(extend.row.dose), 8))
  
  extend.col.dose <- mapply(function(x, y){seq(from = x, to = y,
                                               length.out = len + 2)},
                            col.dose[-nc], col.dose[-1])
  extend.col.dose <- unique(round(c(extend.col.dose), 8))
  
  rownames(extended.scores) <- extend.row.dose
  colnames(extended.scores) <- extend.col.dose
  
  return(extended.scores)
}


# set wd
# I unpacked Sim’s biolog file here, you can do it wherever you want in your computers
setwd("D:/MRC_Postdoc/Pangenomic/biolog/Biolog_metf_40_60")

### Load data ####

data = read_csv('biologdata/Output/Summary.csv', quote = "\"") %>%
  rename(AUC_raw = `750nm_f_AUC`, 
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
         SampleID = paste(Sample, Replicate, Metformin_mM, sep = '_'),
         Sample = as.factor(Sample),
         Strain = as.factor(Strain),
         DrugConc = str_sub(MetaboliteU, -1)) %>% #Change Type column coding
  unite(Type_met, Type, Metformin_mM, remove = FALSE) %>%
  select(SampleID, Metformin_mM, Sample, Strain, Type, Type_met, Replicate, Replicate, Index, Plate, 
         Well, Row, Col, Metabolite, MetaboliteU, EcoCycID, KEGG_ID, Group, Class, AUC_raw) %>%
  mutate_at(c('SampleID', 'Metformin_mM', 'Type_met', 'Metabolite', 'MetaboliteU'), as.factor)


# how many reps, plates, conditions and everything
data %>% group_by(Metformin_mM, Plate, Type) %>% summarise(N = n())

unique(data$Plate)

unique(data$Metformin_mM)

unique(data$Metabolite)



### Summary stats ####

# remove pm10 as it doesn't have MetaboliteU names for some of the compounds
data.sum = data %>% 
  filter(Plate != 'PM10') %>%
  group_by(Metformin_mM, Strain, Sample, Type, Index, Plate, MetaboliteU, Metabolite) %>%
  summarise(Mean = mean(AUC_raw, na.rm = TRUE),
            SD = sd(AUC_raw, na.rm = TRUE))



### Heatmap ####

# lets try and create a heatmap with z-scores by metabolite and plate

met.wide = data %>%
  filter(Plate != 'PM10') %>%
  filter(Metformin_mM %in% c(0,20,40,60)) %>%
  unite(MetID, MetaboliteU, Plate) %>%
  select(SampleID, Metformin_mM, MetID, AUC_raw) %>%
  pivot_wider(names_from = c(SampleID, Metformin_mM), values_from = AUC_raw)

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
             file = here('summary', 'Heatmap.pdf'),
             width = 10, height = 15, useDingbats = FALSE)





# # # # # # # # # # # # # # # #  
### Multi-univariate stats ####
# # # # # # # # # # # # # # # # 




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
  filter(Plate != 'PM10') %>%
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

# similar treatment to get 5-Fluorouracil as controls
drug_ctr = data %>% filter(Metabolite %in% c('5-Fluorouracil')) %>%
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



### increasing the number of compounds ###

drugs = c('Amikacin', 'D-Cycloserine', 'Erythromycin', 'Atropine', 'Orphenadrine',
          'Neomycin', 'Tetracycline', 'Spectinomycin', 'Chloramphenicol', 'Hygromycin B',
          'Geneticin disulfate (G418)', 'Triclosan', 'Procaine', 'Diamide', 
          'Menadione, sodium bisulfite', 'Sodium Dichromate', 'Sodium Tungstate',
          'Sodium Nitrite', 'Alexidine', 'Guanidine hydrochloride', 'Chlorhexidine diacetate',
          'Aztreonam', 'Chlorambucil', 'Sulfamethazine', 'Gallic acid', 'Ethionamide', 
          'Lithium chloride', 'Harmane', 'Benserazide', 'Pridinol', 'Sodium metasilicate',
          'Ethylene Glycol-bis(b-Aminoethyl ether)-NNN`N`-Tetraacetic Acid', 
          '5-Fluorouracil', '5-Fluoroorotic acid', 'Sodium Cyanate', 'Sodium Azide',
          'Methyl viologen')

other_ctr = c('X-b-D-Glucuronide', 'X-Caprylate')
pH = c('pH 5.5', 'pH 6', 'pH 6.5', 'pH 7')


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

# similar treatment to get 5-Fluorouracil as controls
drug_ctr = data %>% filter(Metabolite %in% c('5-Fluorouracil')) %>%
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

num = 1
for (drug in drugs){
  complete.test[complete.test$Drug2 == drug,]$PairIndex = num
  num = num + 1
}



write.csv(complete.test, 'SynergyFinder_test.csv', row.names = FALSE)

complete.test2 = complete.test %>% mutate(Conc1 = ifelse(Conc1 == 40, 4, ifelse(Conc1 == 60, 6, 0)))
write.csv(complete.test2, 'SynergyFinder_test_conc.csv', row.names = FALSE)


# # # # # # # # # # # # # # # # # # # # # 
### SynergyFinder with ALL compounds ####
# # # # # # # # # # # # # # # # # # # # # 

data2 = data %>% 
  filter(Plate !=  'PM10') %>%
  unite(MetPlate, Metabolite, Plate)

drugs = data2 %>% select(MetPlate) %>% t %>% as.vector %>% unique

other_ctr = c('X-b-D-Glucuronide', 'X-Caprylate')
pH = c('pH 5.5', 'pH 6', 'pH 6.5', 'pH 7')


# let's use pH7 as control for drug == 0, metf == (0,40,60)

# substracting drugs and renaming everything
test = data2 %>% filter(MetPlate %in% drugs) %>%
  mutate(DrugConc = str_sub(MetaboliteU, -1),
         Drug1 = 'Metformin',
         ConcUnit = 'A.U.',
         PairIndex = 1) %>%
  select(PairIndex, Metformin_mM, MetPlate, DrugConc, AUC_raw, Drug1, ConcUnit) %>%
  rename(Conc1 = Metformin_mM,
         Conc2 = DrugConc, 
         Drug2 = MetPlate,
         Response = AUC_raw) %>%
  arrange(Drug2, Conc1, Conc2) %>%
  mutate(Conc2 = as.numeric(Conc2)) %>%
  select(PairIndex, Drug1, Drug2, Conc1, Conc2, Response, ConcUnit)

# similar treatment to get 5-Fluorouracil as controls
drug_ctr = data %>% filter(MetaboliteU %in% c('5-Fluorouracil|1')) %>%
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

num = 1
for (drug in drugs){
  complete.test[complete.test$Drug2 == drug,]$PairIndex = num
  num = num + 1
}


# write half of drugs in one file, and other half in a new file
# the web server was crashing with all drugs at once
write.csv(complete.test %>%
            filter(PairIndex < 60), here('summary','SynergyFinder_dataset1.csv'), row.names = FALSE)
write.csv(complete.test %>%
            filter(PairIndex >= 60 & PairIndex < 120), here('summary','SynergyFinder_dataset2.csv'), row.names = FALSE)
write.csv(complete.test %>%
            filter(PairIndex >= 120 & PairIndex < 180), here('summary','SynergyFinder_dataset3.csv'), row.names = FALSE)
write.csv(complete.test %>%
            filter(PairIndex >= 180), here('summary','SynergyFinder_dataset4.csv'), row.names = FALSE)



# # # # # # # # # # # # #
# WORKING ZONE: DANGER #
# # # # # # # # # # # # #


### SynergyFinder R version ####
library(synergyfinder)
library(plotly)

data("mathews_screening_data")
head(mathews_screening_data)

# reshape and pre-process the input data to a 
# dose-response matrix format for further analysis:

set.seed(1)
dose.response.mat <- ReshapeData(mathews_screening_data,
                                 data.type = "viability",
                                 impute = TRUE,
                                 noise = TRUE,
                                 correction = "non")


# plot heatmaps
PlotDoseResponse(dose.response.mat)

# calculate synergy score
synergy.score <- CalculateSynergy(data = dose.response.mat,
                                  method = "ZIP")


PlotSynergy(synergy.score, type = "all")

data("mathews_screening_data")
data <- ReshapeData(mathews_screening_data)
response.mat <- data$dose.response.mats[[1]]
ZIP.score <- ZIP(response.mat)

rwn = rownames(ZIP.score)

cosa = as_tibble(ZIP.score) %>%
  mutate(drug.conc1 = nam) %>%
  pivot_longer(cols = `0`:`50`, names_to = 'drug.conc2') %>%
  mutate_at(c('drug.conc1', 'drug.conc2'), as.factor)
  
  

fig <- plot_ly(data = cosa,x = ~drug.conc1, y = ~drug.conc2, z = ~value)
fig <- fig %>% add_surface()

fig




fig <- plot_ly(z = ~ZIP.score)
fig <- fig %>% add_surface()

fig



ZIP2 = ExtendedScores(ZIP.score, len = 3)

# fig <- plot_ly(z = ~ZIP2)
# fig <- fig %>% add_surface()
# 
# fig

# profesional plot

x.conc = colnames(ZIP.score)
y.conc = rownames(ZIP.score)

subX = seq.int(min(x.conc), max(x.conc), length.out = nrow(ZIP2))
subY = seq.int(max(y.conc), min(y.conc), length.out = ncol(ZIP2))

plot_ly(z = rotate(rotate(rotate(ZIP2))), x = subX, y = subY, type = "surface", 
            contours = list( y = list(show = !0, width = 1, highlightwidth = 2, usecolormap = F),  
                             x = list(show = !0, width = 1, highlightwidth = 2, usecolormap = F)), 
            hoverinfo = "z+name", name = "d - score", 
            colorbar = list(outlinewidth = 0, title = "\U03B4 - score", len = 0.24, thickness = 19, xpad = 3, showticklabels = !0, 
                            titlefont = 9, outlinewidth = 0.3, tickcolor = "#fff", tickfont = list(size = 9), ticks = "inside"), 
            colorscale = list(c(0, "rgb(0, 247, 0)"), c(0.5, "rgb(247, 247, 247)"), c(1, "rgb(247, 0, 0)")), 
            cauto = F, cmin = -35, cmax = 35, 
            # contour = list(show = !1, color = "#222")
        )

PlotSynergy(synergy.score, type = "all")

mean(ZIP2)



########################

data("mathews_screening_data")
data <- ReshapeData(mathews_screening_data)
response.mat <- data$dose.response.mats[[1]]
ZIP.score <- ZIP(response.mat)



# lets try to reshape one of the drugs from my data

syn.test = complete.test %>%
  filter(Drug2 == 'Orphenadrine') %>%
  rename(blockID = PairIndex,
         drug_col = Drug1,
         drug_row = Drug2,
         conc_c = Conc1,
         conc_r = Conc2,
         response = Response,
         conc_c_unit = ConcUnit) %>%
  mutate(conc_r_unit = conc_c_unit)

data <- ReshapeData(syn.test)


response.mat <- data$adjusted.response.mats[[1]]
ZIP.score <- ZIP(response.mat)


# let's use pH7 as control for drug == 0, metf == (0,40,60)

# substracting drugs and renaming everything
test = data %>% filter(Metabolite %in% drugs) %>%
  mutate(DrugConc = str_sub(MetaboliteU, -1),
         Drug1 = 'Metformin',
         ConcUnit = 'A.U.',
         PairIndex = 1) %>%
  select(PairIndex, Metformin_mM, Metabolite, Replicate, DrugConc, AUC_raw, Drug1, ConcUnit) %>%
  rename(Conc1 = Metformin_mM,
         Conc2 = DrugConc, 
         Drug2 = Metabolite,
         Response = AUC_raw) %>%
  arrange(Drug2, Conc1, Conc2) %>%
  mutate(Conc2 = as.numeric(Conc2)) %>%
  select(PairIndex, Drug1, Drug2, Replicate, Conc1, Conc2, Response, ConcUnit)

# similar treatment to get 5-Fluorouracil as controls
drug_ctr = data %>% filter(Metabolite %in% c('5-Fluorouracil')) %>%
  mutate(DrugConc = str_sub(MetaboliteU, -1),
         DrugConc = 0,
         Drug1 = 'Metformin',
         ConcUnit = 'A.U.',
         PairIndex = 1) %>%
  select(PairIndex, Metformin_mM, Metabolite, Replicate, DrugConc, AUC_raw, Drug1, ConcUnit) %>%
  rename(Conc1 = Metformin_mM,
         Conc2 = DrugConc,
         Drug2 = Metabolite,
         Response = AUC_raw) %>%
  arrange(Drug2, Conc1, Conc2) %>%
  select(PairIndex, Drug1, Conc1, Conc2, Replicate, Response, ConcUnit) 

# correct for number showing in replicate column
for (conc in c(0,20,40,60)){
  drug_ctr[drug_ctr$Conc1 == conc,]$Replicate = 1:length(drug_ctr[drug_ctr$Conc1 == 60,]$Replicate)
}

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

num = 1
for (drug in drugs){
  complete.test[complete.test$Drug2 == drug,]$PairIndex = num
  num = num + 1
}


complete.test


syn.test = complete.test %>%
  filter(Drug2 == 'Orphenadrine') %>%
  rename(blockID = PairIndex,
         drug_col = Drug1,
         drug_row = Drug2,
         conc_c = Conc1,
         conc_r = Conc2,
         response = Response,
         conc_c_unit = ConcUnit,
         replicate = Replicate) %>%
  mutate(conc_r_unit = conc_c_unit)

data <- ReshapeData(syn.test)


response.mat <- data$adjusted.response.mats[[1]]
ZIP.score <- ZIP(response.mat)



# # # # # # # # # # # #
### Synergy plots ####
# # # # # # # # # # # #

### Bliss
# read excel files

bliss_sc_files = c('result_Bliss_scores_1.xlsx', 'result_Bliss_scores_2.xlsx', 
                   'result_Bliss_scores_3.xlsx', 'result_Bliss_scores_4.xlsx')
bliss.sc = tibble()
for (file in bliss_sc_files) {
  temp = read_excel(paste0('summary/SynergyFinder/', file))
  bliss.sc = rbind(bliss.sc, temp)
}

# remove tail in name, and reorder by most synergistic
bliss.sc = bliss.sc %>% rename(CI = `95% CI`) %>%
  arrange(Synergy.score) %>%
  mutate(Drug.combination = substr(Drug.combination, 1, nchar(Drug.combination) - 12),
         Drug.combination = fct_reorder(Drug.combination, desc(Synergy.score))) # remove Metformin from names

bliss.sc %>%
  ggplot(aes(x = Drug.combination, y = Synergy.score)) +
  geom_hline(yintercept = 0, size = 1, colour = 'grey50') +
  geom_point() +
  geom_errorbar(aes(x = Drug.combination, ymin = Synergy.score - CI, ymax = Synergy.score + CI)) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust  = 1)) 


# save plot
dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'Bliss_ordered.pdf'),
             width = 34, height = 10, useDingbats = FALSE)

write.csv(bliss.sc, here('summary','Bliss_scores_corrected.csv'), row.names = F)


### ZIP


zip_sc_files = c('result_ZIP_scores_1.xlsx', 'result_ZIP_scores_3.xlsx', 
                   'result_ZIP_scores_4.xlsx')
zip.sc = tibble()
for (file in zip_sc_files) {
  temp = read_excel(paste0('summary/SynergyFinder/', file))
  zip.sc = rbind(zip.sc, temp)
}


# remove tail in name, and reorder by most synergistic
zip.sc = zip.sc %>% rename(CI = `95% CI`) %>%
  arrange(Synergy.score) %>%
  mutate(Drug.combination = substr(Drug.combination, 1, nchar(Drug.combination) - 12),
         Drug.combination = fct_reorder(Drug.combination, desc(Synergy.score))) # remove Metformin from names

zip.sc %>%
  ggplot(aes(x = Drug.combination, y = Synergy.score)) +
  geom_hline(yintercept = 0, size = 1, colour = 'grey50') +
  geom_point() +
  geom_errorbar(aes(x = Drug.combination, ymin = Synergy.score - CI, ymax = Synergy.score + CI)) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust  = 1)) 


# save plot
dev.copy2pdf(device = cairo_pdf,
             file = here('summary', 'ZIP_ordered_INCOMPLETE.pdf'),
             width = 28, height = 10, useDingbats = FALSE)

