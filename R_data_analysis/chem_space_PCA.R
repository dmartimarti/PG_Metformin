# libraries

library(tidyverse) # master library to deal with data frames
library(readxl) # read xlsx or xls files
library(ggrepel) # ggplot add-on, to plot names that don't collapse in same position
library(FactoMineR) # for PCA
library(factoextra) # for PCA
library(here) # usefull to save plots in folders given a root
library(viridis) # color palette package
# library(ComplexHeatmap) # yeah, complex heatmaps
library(plotly)


# load necessary files to plot PCA

# biolog metabolites
biolog = read_csv('Biolog_metabolites.csv')

# drugbank coordinates
drugbank.coord = read.delim('.\\Data_drug_bacteria_gene_mapping\\Input\\drugbank_pca.txt', header = FALSE)
drugbank.coord2 = read.delim('.\\Data_drug_bacteria_gene_mapping\\Input\\pharmacon_pca.txt', header = FALSE)
drugbank = read_csv('.\\Data_drug_bacteria_gene_mapping\\Input\\drugbank_approved_MW_150_1000_functional_groups_all.csv') %>% select(-Number_of_carboxylic_acids_1)

plates = c("PM11C", "PM12B", "PM13B", "PM14A", "PM15B", "PM16A", "PM17A", "PM18C", "PM19" , "PM20B")
biolog.drug = biolog %>% filter(Plate %in% plates)
drugs = biolog.drug %>% select(Metabolite) %>% unique %>% t %>% as.character

length(drugs[toupper(drugs) %in% drugbank.coord2$V1])

drugs[!drugs %in% unique(drugbank$GENERIC_NAME)]

# so many compounds without match between databases, let's change the giant csv file

drugbank.coord['drugs'] = drugbank$GENERIC_NAME

test = drugbank.coord
rownames(test) = drugbank$GENERIC_NAME

fig <- plot_ly(data = test, x = ~V1, y = ~V2, z = ~V3, text = rownames(test),
               marker = list(size = 4,
                             color = '#808080'))
fig

##

drugbank.coord['Drugbank'] = drugbank$GENERIC_NAME

test2 = drugbank.coord

test2 = test2 %>%
  mutate(biolog = ifelse(Drugbank %in% drugs, Drugbank, NA),
         Category = ifelse(Drugbank %in% drugs, 'Biolog', 'Drugbank'),
         Category = as.factor(Category)) 

rownames(test2) = test2$Drugbank

fig <- plot_ly(data = test2, x = ~V1, y = ~V2, z = ~V3, text = rownames(test2), 
               color = ~Category, colors = c('#FF0900', '#B8B8B8'),
               marker = list(size = 3))
fig




### let's get the missing drugs (as many as we can, at least)
missing = drugs[!drugs %in% unique(drugbank$GENERIC_NAME)]

# first, separate KEGG ids to generate a new column to the biolog data frame with pubchem ids
kegg_ids = biolog %>%
  # filter(Metabolite %in% missing) %>% 
  select(KEGG_ID) %>% t %>% as.character %>% unique 
# remove NAs
kegg_ids = kegg_ids[complete.cases(kegg_ids)]
# write list of genes
write.table(kegg_ids, 'KEGG_IDs_biolog.txt', quote = F, col.names = F, row.names = F)
## go here to convert: http://csbg.cnb.csic.es/mbrole/conversion.jsp

# read generated list and merge with biolog
kegg2pub =  read.csv("D:/MRC_Postdoc/Pangenomic/Chem_space/KEGG2PubChemIDs.txt")
biolog = biolog %>% left_join(kegg2pub)


drugs_missing = biolog %>%
  filter(Metabolite %in% missing) %>%
  select(PubChem_ID) %>% t %>% as.character %>% unique 
# remove NAs
drugs_missing = drugs_missing[complete.cases(drugs_missing)]
write.table(drugs_missing, 'PubChem_ID_missing_metabolites.txt', quote = F, col.names = F, row.names = F)

# after downloading them from pubchem, the file is named as: Metab_structures_missing_biolog.sdf

# file with chem fingerprints from pubchem compounds (~100 more, not bad)
biolog_metabolites_PubChem = read_csv("biolog_metabolites_PubChem.csv") %>% select(-X1, -index, -Number_of_carboxylic_acids_1, -ID) %>%
  rename(PubChem_ID = PUBCHEM_COMPOUND_CID) 


dummy = biolog %>%
  filter(Metabolite %in% missing) %>%
  select(Metabolite, PubChem_ID) %>% unique

biolog_metabolites_PubChem = biolog_metabolites_PubChem %>% left_join(dummy, by = 'PubChem_ID') %>%
  select(PubChem_ID, Metabolite, everything()) %>% rename(GENERIC_NAME = Metabolite)


# bind rows
complete = bind_rows(drugbank %>% 
                       mutate(Category = 'Drugbank'), 
                     biolog_metabolites_PubChem) %>% 
  distinct(GENERIC_NAME, .keep_all = TRUE)

# save data table into csv
write.csv(complete, 'Complete_chem_fingerprints.csv')



mat = complete %>% 
  filter(GENERIC_NAME != 'Tannic acid') %>%
  # select(-Number_of_carboxylic_acids_1) %>%
  select(Number_of_aliphatic_carboxylic_acids:Number_of_urea_groups) %>% as.matrix
rownames(mat) = complete %>% filter(GENERIC_NAME != 'Tannic acid') %>% select(GENERIC_NAME) %>% t %>% as.character


res.pca = PCA((mat), scale.unit = TRUE, ncp = 5, graph = F)
ind = get_pca_ind(res.pca)
ind_df = data.frame(ind$coord[,1], ind$coord[,2], ind$coord[,3])
colnames(ind_df) = c('Dim1', 'Dim2', 'Dim3')


## plot with names 
test2 = ind_df %>% tibble 
test2['Drugbank'] = complete %>% filter(GENERIC_NAME != 'Tannic acid') %>% select(GENERIC_NAME) %>% t %>% as.character

test2 = test2 %>%
  mutate(biolog = ifelse(Drugbank %in% drugs, Drugbank, NA),
         Category = ifelse(Drugbank %in% drugs, 'Biolog', 'Drugbank'),
         Category = as.factor(Category)) 


fig = plot_ly(data = test2,  x = ~Dim1, y = ~Dim2, z = ~Dim3, text = test2$Drugbank, 
               color = ~Category, colors = c('#FF0900', '#B8B8B8'),
               marker = list(size = 3))
fig


### let's try t-SNE

library(Rtsne)

# scale the data first prior to running t-SNE

mat.tsne = complete %>% 
  # filter(GENERIC_NAME != 'Tannic acid') %>%
  # select(-Number_of_carboxylic_acids_1) %>%
  select(Number_of_aliphatic_carboxylic_acids:Number_of_urea_groups) %>% as.matrix
rownames(mat.tsne) = complete %>% 
  # filter(GENERIC_NAME != 'Tannic acid') %>% 
  select(GENERIC_NAME) %>% t %>% as.character

tsne = Rtsne(mat.tsne, check_duplicates = FALSE, pca = FALSE, num_threads = 12,
             normalize = FALSE, max_iter = 2000,
             perplexity = 40, theta = 0.2, dims = 3) 

# generate data frame from tnse results
tsne.df = data.frame(tsne$Y)
colnames(tsne.df) = c('Dim1', 'Dim2', 'Dim3')

tsne.df['Drugbank'] = complete %>% 
  # filter(GENERIC_NAME != 'Tannic acid') %>% 
  select(GENERIC_NAME) %>% t %>% as.character

tsne.df = tsne.df %>%
  mutate(biolog = ifelse(Drugbank %in% drugs, Drugbank, NA),
         Category = ifelse(Drugbank %in% drugs, 'Biolog', 'Drugbank'),
         Opacity = ifelse(Category == 'Biolog', 1, 0.1),
         Opacity = as.numeric(Opacity),
         Category = as.factor(Category)) 


fig = plot_ly(data = tsne.df,  x = ~Dim1, y = ~Dim2, z = ~Dim3, text = tsne.df$Drugbank, 
               color = ~Category, colors = c('#FF0900', '#B8B8B8'),
               marker = list(size = 4))
fig


# let's try to modify the opacity of grey points
fig = plot_ly() 

fig = fig %>%
  add_trace(
    data = tsne.df %>% filter(Category == 'Drugbank'),
    name = 'Drugbank',
    x = ~Dim1, y = ~Dim2, z = ~Dim3, text = tsne.df[tsne.df$Category == 'Drugbank',]$Drugbank, 
    marker = list(
      color = '#B8B8B8', 
      opacity = 0.4, # OPACITY
      size = 4)
    )

fig = fig %>%
  add_trace(
    data = tsne.df %>% filter(Category == 'Biolog'),
    name = 'Biolog',
    x = ~Dim1, y = ~Dim2, z = ~Dim3, text = tsne.df[tsne.df$Category == 'Biolog',]$Drugbank, 
    marker = list(
      color = '#FF0900', 
      size = 4)
  )
 
fig


### and now a version with a slider for the opacity

fig = plot_ly(type = 'scatter3d', mode = 'markers') %>% 
  layout(sliders = list(
    list(
      active = 0, 
      currentvalue = list(prefix = "Opacity: "), 
      pad = list(t = 60), 
      steps = steps))
    )

# create steps for slider
steps = list(
  list(args = list("marker.opacity", 0.3), 
       label = "0.3", 
       method = "restyle", 
       value = "1"
  ),
  list(args = list("marker.opacity", 0.6), 
       label = "0.6", 
       method = "restyle", 
       value = "2"
  ),
  list(args = list("marker.opacity", 1), 
       label = "1", 
       method = "restyle", 
       value = "3"
  )
)

fig = fig %>%
  add_trace(
    data = tsne.df %>% filter(Category == 'Drugbank'),
    name = 'Drugbank',
    x = ~Dim1, y = ~Dim2, z = ~Dim3, text = tsne.df[tsne.df$Category == 'Drugbank',]$Drugbank, 
    marker = list(
      color = '#B8B8B8', 
      opacity = 0.4, # OPACITY
      size = 4)
  ) 


fig = fig %>%
  add_trace(
    data = tsne.df %>% filter(Category == 'Biolog'),
    name = 'Biolog',
    x = ~Dim1, y = ~Dim2, z = ~Dim3, text = tsne.df[tsne.df$Category == 'Biolog',]$Drugbank, 
    marker = list(
      color = '#FF0900', 
      size = 4)
  )

fig 


### Completing biolog drugs ####


biolog.table = biolog %>% 
  filter(Plate %in% plates) %>% 
  select(Metabolite, Plate, KEGG_ID, CAS_ID, PubChem_ID) %>% 
  unique

# write table
write.csv(biolog.table, 'biolog_pubchem_IDs.csv')

# I LOOKED FOR THE MISSING COMPOUNDS AND ADDED THEM
# ALSO, I ADDED A NEW COLUMN WITH PUBCHEM SUBSTANCE
# AS 3 COMPOUNDS WERE NOT FOUND IN PUBCHEM

# load AGAIN the file, it has been updated
biolog_pubchem_IDs = read_csv("biolog_pubchem_IDs.csv")

# store the list of pubchem IDs
pubchem_IDs = biolog_pubchem_IDs %>% select(PubChem_ID) %>% t %>% as.character 
pubchem_IDs = pubchem_IDs[complete.cases(pubchem_IDs)]

# store the list of pubchem substances
pubchem_subs = biolog_pubchem_IDs %>% select(PubChem_Substance) %>% t %>% as.character 
pubchem_subs = pubchem_subs[complete.cases(pubchem_subs)]

# write pubchem ids
write.table(pubchem_IDs, 'PubChem_ID_missing_metabolites.txt', quote = F, col.names = F, row.names = F)
# write pubchem substances
write.table(pubchem_subs, 'PubChem_Subs_missing_metabolites.txt', quote = F, col.names = F, row.names = F)

# AFTER THAT YOU MUST RUN THE PYTHON SCRIPT FOR EVERY METABOLITE
# NOW IT IS ONLY WORKING THE PUBCHEM ID, SUBSTANCES ARE NOT ACCEPTED BY THE SCRIPT, NEED TO CHANGE THAT


# file with chem fingerprints from pubchem compounds (~100 more, not bad)
biolog_metabolites_PubChem = read_csv("biolog_metabolites_PubChem.csv") %>% select(-X1, -index, -Number_of_carboxylic_acids_1, -ID) %>%
  rename(PubChem_ID = PUBCHEM_COMPOUND_CID) 

dummy = biolog_pubchem_IDs %>%
  select(Metabolite, PubChem_ID) %>% unique

biolog_metabolites_PubChem = biolog_metabolites_PubChem %>% left_join(dummy, by = 'PubChem_ID') %>%
  select(PubChem_ID, Metabolite, everything()) %>% rename(GENERIC_NAME = Metabolite)

# add metformin to the database
biolog_metabolites_PubChem[biolog_metabolites_PubChem$PubChem_ID == 4091,]$GENERIC_NAME <- 'Metformin'

# add a new category for biolog compounds
biolog_metabolites_PubChem = biolog_metabolites_PubChem %>% 
  mutate(Category = 'Biolog')

# bind rows
complete = bind_rows(drugbank %>% 
                       mutate(Category = 'Drugbank'), 
                     biolog_metabolites_PubChem) %>% 
  distinct(GENERIC_NAME, .keep_all = TRUE)


# save data table into csv
write.csv(complete, 'Complete_chem_fingerprints.csv')


## t-SNE Complete (Almost) ####

# bind rows
complete = bind_rows(drugbank %>%
                       mutate(Category = 'Drugbank'),
                     biolog_metabolites_PubChem) %>%
  distinct(GENERIC_NAME, .keep_all = TRUE)

library(Rtsne)

# scale the data first prior to running t-SNE

mat.tsne = complete %>% 
  # filter(GENERIC_NAME != 'Tannic acid') %>%
  # select(-Number_of_carboxylic_acids_1) %>%
  select(Number_of_aliphatic_carboxylic_acids:Number_of_urea_groups) %>% as.matrix
rownames(mat.tsne) = complete %>% 
  # filter(GENERIC_NAME != 'Tannic acid') %>% 
  select(GENERIC_NAME) %>% t %>% as.character

tsne = Rtsne(mat.tsne, check_duplicates = FALSE, pca = FALSE, num_threads = 12,
             normalize = FALSE, max_iter = 2000,
             perplexity = 40, theta = 0.2, dims = 3) 

# generate data frame from tnse results
tsne.df = data.frame(tsne$Y)
colnames(tsne.df) = c('Dim1', 'Dim2', 'Dim3')

tsne.df['Drugbank'] = complete %>% 
  # filter(GENERIC_NAME != 'Tannic acid') %>% 
  select(GENERIC_NAME) %>% t %>% as.character


tsne.df = tsne.df %>%
  mutate(biolog = ifelse(Drugbank %in% drugs, Drugbank, NA),
         Category = ifelse(Drugbank %in% drugs, 'Biolog', 'Drugbank'),
         Category = ifelse(Drugbank == 'Metformin', 'Metformin', Category),
         Category = as.factor(Category)) 


# let's try to modify the opacity of grey points
fig = plot_ly() 

fig = fig %>%
  add_trace(
    data = tsne.df %>% filter(Category == 'Drugbank'),
    name = 'Drugbank',
    x = ~Dim1, y = ~Dim2, z = ~Dim3, text = tsne.df[tsne.df$Category == 'Drugbank',]$Drugbank, 
    marker = list(
      color = '#B8B8B8', 
      opacity = 0.4, # OPACITY
      size = 4)
  )

fig = fig %>%
  add_trace(
    data = tsne.df %>% filter(Category == 'Biolog'),
    name = 'Biolog',
    x = ~Dim1, y = ~Dim2, z = ~Dim3, text = tsne.df[tsne.df$Category == 'Biolog',]$Drugbank, 
    marker = list(
      color = '#FF0900', 
      size = 4)
  )

fig = fig %>%
  add_trace(
    data = tsne.df %>% filter(Category == 'Metformin'),
    name = 'Metformin',
    x = ~Dim1, y = ~Dim2, z = ~Dim3, text = tsne.df[tsne.df$Category == 'Metformin',]$Drugbank, 
    marker = list(
      color = '#3E4AFF', 
      size = 6)
  )


fig

# # # # # # # # # # # # # # # # # #
### Synergy/antagonistic plot ####
# # # # # # # # # # # # # # # # # #


bliss = read_csv('D:/MRC_Postdoc/Pangenomic/biolog/Biolog_metf_40_60/summary/Bliss_scores_corrected.csv')

bliss = bliss %>% separate(Drug.combination, c('Drug', 'Plate'), sep = '_')

# initiate the threshold
thrs = 4

### BE CAREFUL AND CHOSE A SPECIFIC WAY TO SEPARATE COMPOUNDS
# classify by CI_low
bliss.tsne = bliss %>% 
  mutate(CI_low = abs(Synergy.score) - abs(CI),
         CI_up = abs(Synergy.score) + abs(CI),
         Direction = ifelse(CI_low > thrs & Synergy.score < 0, 'Antagonistic', 
                            ifelse(CI_low > thrs & Synergy.score > 0, 'Synergistic', 'Neutral')),
         Drug = fct_reorder(Drug, desc(Synergy.score))) %>%
  rename(Drugbank = Drug)



tsne.df = tsne.df %>% left_join(bliss.tsne)




# let's try to modify the opacity of grey points
fig = plot_ly() 

fig = fig %>%
  add_trace(
    data = tsne.df %>% filter(Category == 'Drugbank'),
    name = 'Drugbank',
    x = ~Dim1, y = ~Dim2, z = ~Dim3, text = (tsne.df %>% filter(Category == 'Drugbank'))$Drugbank, 
    marker = list(
      color = '#B8B8B8', 
      opacity = 0.4, # OPACITY
      size = 4)
  )

fig = fig %>%
  add_trace(
    data = tsne.df %>% filter(Direction == 'Antagonistic'),
    name = 'Antagonistic',
    x = ~Dim1, y = ~Dim2, z = ~Dim3, text = (tsne.df %>% filter(Direction == 'Antagonistic'))$Drugbank ,
    marker = list(
      color = '#FF0900', 
      size = 6)
  )

fig = fig %>%
  add_trace(
    data = tsne.df %>% filter(Direction == 'Synergistic'),
    name = 'Synergistic',
    x = ~Dim1, y = ~Dim2, z = ~Dim3, text = (tsne.df %>% filter(Direction == 'Synergistic'))$Drugbank, 
    marker = list(
      color = '#057D33', 
      size = 6)
  )

fig = fig %>%
  add_trace(
    data = tsne.df %>% filter(Direction == 'Neutral'),
    name = 'Neutral',
    x = ~Dim1, y = ~Dim2, z = ~Dim3, text = (tsne.df %>% filter(Direction == 'Neutral'))$Drugbank, 
    marker = list(
      color = '#4446E0', 
      size = 6,
      opacity = 1)
  )

fig = fig %>%
  add_trace(
    data = tsne.df %>% filter(Category == 'Metformin'),
    name = 'Metformin',
    x = ~Dim1, y = ~Dim2, z = ~Dim3, text = (tsne.df %>% filter(Category == 'Metformin'))$Drugbank, 
    marker = list(
      color = '#000000', 
      size = 6)
  )


fig


### ML Prediction ####

bliss.drugs = bliss.tsne %>% select(Drugbank) %>% t %>% as.vector

# compounds from our biolog
forest_train = complete %>% filter(GENERIC_NAME %in% bliss.drugs) %>%
  select(GENERIC_NAME, Number_of_aliphatic_carboxylic_acids:Number_of_urea_groups) %>%
  left_join(bliss.tsne %>% rename(GENERIC_NAME = Drugbank) %>%
              select(GENERIC_NAME, Direction))

write.table(forest_train, 'ML_training_set.txt', quote = F,  row.names = F, sep = '\t')

# save the rest of the compounds from the big database
forest_predict = complete %>% filter(!GENERIC_NAME %in% bliss.drugs) %>%
  select(GENERIC_NAME, Number_of_aliphatic_carboxylic_acids:Number_of_urea_groups) 

write.csv(forest_predict, 'ML_predict.csv', quote = F,  row.names = F)

## after having trained my model, I have predicted values for the other compounds using KNN algorithm

ML_predicted_KNN = read_excel("ML_predicted_KNN.xlsx")


ML_predicted_KNN = ML_predicted_KNN %>%
  # select(GENERIC_NAME, Direction) %>%
  rename(Drugbank = GENERIC_NAME)

# generate data frame from tnse results
tsne.df = data.frame(tsne$Y)
colnames(tsne.df) = c('Dim1', 'Dim2', 'Dim3')

tsne.df['Drugbank'] = complete %>% 
  # filter(GENERIC_NAME != 'Tannic acid') %>% 
  select(GENERIC_NAME) %>% t %>% as.character


tsne.df = tsne.df %>%
  mutate(biolog = ifelse(Drugbank %in% drugs, Drugbank, NA),
         Category = ifelse(Drugbank %in% drugs, 'Biolog', 'Drugbank'),
         Category = ifelse(Drugbank == 'Metformin', 'Metformin', Category),
         Category = as.factor(Category)) 

tsne.df = tsne.df %>% 
  left_join(bliss.tsne %>%
              bind_rows(ML_predicted_KNN))


### plot figure



# let's try to modify the opacity of grey points
fig = plot_ly() 

fig = fig %>%
  add_trace(
    data = tsne.df %>% filter(Direction == 'Antagonistic'),
    name = 'Antagonistic',
    x = ~Dim1, y = ~Dim2, z = ~Dim3, text = (tsne.df %>% filter(Direction == 'Antagonistic'))$Drugbank ,
    marker = list(
      color = '#FF0900', 
      size = 6)
  )

fig = fig %>%
  add_trace(
    data = tsne.df %>% filter(Direction == 'Synergistic'),
    name = 'Synergistic',
    x = ~Dim1, y = ~Dim2, z = ~Dim3, text = (tsne.df %>% filter(Direction == 'Synergistic'))$Drugbank, 
    marker = list(
      color = '#057D33', 
      size = 6)
  )

fig = fig %>%
  add_trace(
    data = tsne.df %>% filter(Direction == 'Neutral'),
    name = 'Neutral',
    x = ~Dim1, y = ~Dim2, z = ~Dim3, text = (tsne.df %>% filter(Direction == 'Neutral'))$Drugbank, 
    marker = list(
      color = '#4446E0', 
      size = 6,
      opacity = 1)
  )

fig = fig %>%
  add_trace(
    data = tsne.df %>% filter(Category == 'Metformin'),
    name = 'Metformin',
    x = ~Dim1, y = ~Dim2, z = ~Dim3, text = (tsne.df %>% filter(Category == 'Metformin'))$Drugbank, 
    marker = list(
      color = '#000000', 
      size = 6)
  )


fig


tsne.df %>%
  filter(Direction == 'Synergistic') %>%
  drop_na(Antag_prob)


