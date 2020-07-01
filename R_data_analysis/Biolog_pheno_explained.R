# This script is to star re-analising Pov's data. I will try to explain 
# how analyses work in R, what is the code required to analyse data in 
# different ways (PCAs, stats, heatmaps, z-scores, clustering, enrichment...)


# Every script in R usually starts loading all the required libraries
# This can change some times when some of your required libraries
# can have conflicts with others. I will explain this later if it's
# needed. 

# libraries

library(tidyverse) # master library to deal with data frames
library(readxl) # read xlsx or xls files
library(ggrepel) # ggplot add-on, to plot names that don't collapse in same position
library(FactoMineR) # for PCA
library(factoextra) # for PCA
library(here) # usefull to save plots in folders given a root
library(viridis) # color palette package
library(ComplexHeatmap) # yeah, complex heatmaps

# in the case you don't have some of these packages, install them with:

install.packages(c('tidyverse', 'ggrepel', 'FactoMineR', 'factoextra', 'here', 'viridis'))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")

# I also add this option to my scripts, to adapt R window size to my screen size
# change the width according to your screen

options(width = 220)


#####################
### Read the data ###
#####################

# In EVERY script I write, I also write the path of my analysis, data, etc
# This is super important as it will save you time in the future if you 
# come back to the script in some time and dont remember where the f***k is everything :)

# get working directory = getwd()
# set working directory, you guess it = setwd()
path = "D:/MRC_Postdoc/Pangenomic/Biolog_phenotyping"

# let's read the data

data = read_csv('All_data/Summary.csv', quote = "\"") %>%
  rename(AUC_raw = `750nm_f_AUC24`) %>% # `750nm_f_logAUC` data column is what we need for logAUC values
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
         SampleID = paste(Sample, Replicate, sep = '_'),
         Sample = as.factor(Sample),
         Strain = as.factor(Strain)) %>% #Change Type column coding
  select(SampleID, Sample, Strain, Type, Inoculum, Sugar_20mM, Replicate, Uracil_uM, Metformin_mM, Replicate, Index, Plate, 
         Well, Row, Col, Metabolite, MetaboliteU, EcoCycID, KEGG_ID, Group, Class, AUC_raw)

# you can access different segments of data to see all variables we have

# plates 
unique(data$Plate)

# strains
unique(data$Strain)

# sample
unique(data$Sample)

# inoculum
unique(data$Inoculum)

# uracil
unique(data$Uracil_uM)

# it seems that we have several plates in this data set, and not every one should be treated 
# equal: for example, PM10 to PM20 don't seem to have any 'Negative Control', se we need to split
# our data set into two different datasets

# data belonging to PM1 to PM4, with metabolites
met.data = data %>% filter(Plate %in% c('PM1', 'PM2A', 'PM3B', 'PM4A'), Inoculum == 'LB')

# check again variables in the data frame
unique(met.data$Inoculum)

# check uracil concentrations used in LB media experiments
unique(met.data %>% filter(Inoculum == 'LB') %>% select(Uracil_uM))

# check uracil concentrations used in LB media experiments
unique(met.data %>% filter(Inoculum == 'LB') %>% select(Sugar_20mM))
# several sugar sources, but...

# if we filter by PM1 and PM2, we can see that no additional sugar source was added. 
unique(met.data %>% 
         filter(Plate %in% c('PM1', 'PM2A')) %>% 
         select(Sugar_20mM)
       )


# lets replace met.data with only PM1 and PM2, and lets continue with only those plates
# also, lets remove the negative control values to the other values, to normalise them

met.data = data %>% 
  filter(Plate %in% c('PM1', 'PM2A'), Inoculum == 'LB') %>%
  group_by(Strain, Type, Plate, Replicate, Group, Inoculum, Sugar_20mM, Uracil_uM, Metformin_mM) %>%
  mutate(AUC = AUC_raw / AUC_raw[Metabolite == 'Negative Control'],
         logAUC = log2(AUC)) %>% 
  ungroup


# generate summary statistics for each metabolite in each condition
met.data.sum = met.data %>% 
  select(-Sugar_20mM, -Inoculum) %>%
  group_by(Strain, Type, Plate, Group, Uracil_uM, Metformin_mM, Metabolite, MetaboliteU) %>%
  summarise(Mean = mean(logAUC, na.rm = TRUE),
            SD = sd(logAUC, na.rm = TRUE)) %>%
  ungroup



# so, wait a moment to see what we have done
# For now, we loaded the data, took a look to the variables, and decided to split the set into
# a smaller dataset containing only the PM1 and PM2 plates (PM3 and 4 have other C sources)
# Also, to normalise ODs against neg control, we have divided every value by the neg control,
# and then, after that, calulated the log2 of that value. This will give us the idea of Fold Change (FC)
# I havent done that before because there are some neg controls that are small, and the log2 of a small
# number (<0) is negative. Thus, the normalisation would have worded in the wrong way


##########################
### Heatmap of z-score ###
##########################

# lets try and create a heatmap with z-scores by metabolite and plate
# this will imply some data shaping and manipulation

# first step: make data wide. Select uracil == 40 and remove Neg control
# select only the variables you want to keep, and make it wide
met.wide = met.data %>%
  filter(Uracil_uM == 40) %>%
  filter(Metabolite != 'Negative Control') %>%
  select(SampleID, MetaboliteU, logAUC) %>%
  pivot_wider(names_from = c(SampleID), values_from = logAUC)

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






















