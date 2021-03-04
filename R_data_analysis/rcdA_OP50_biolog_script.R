
# libraries ---------------------------------------------------------------


library(tidyverse)
library(readr)
library(broom)
library(gtools)
library(openxlsx)

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
test = data %>% filter(MetaboliteU == 'Amikacin|2')

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

# this piece of code works in a very simple way:
# take the dataset, group by unique metabolite, nest the data into blocks, and 
# calculate the t-test per group
# unnest everything, tidy everything a bit, and just enjoy your new stats :)


stats = data %>% 
  group_by(MetaboliteU) %>%                               # group by metabolites
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

# look that everything is ok
stats

# save everything into an excel file
list_of_tables = list(
  stats_R = stats
)

write.xlsx(list_of_tables, 'Multi_univariate_stats.xlsx')




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


ggsave('pca_classes.pdf', height = 8,width = 9)


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


ggsave('umap_classes.pdf', height = 8,width = 9)





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

met.wide = data.sum %>%
  ungroup %>% 
  mutate(DrugConc = str_sub(MetaboliteU, -1),
         Plate = str_sub(Index, 1,4)) %>% 
  unite(MetID, Metabolite, Plate) %>%
  unite(ID, Strain, Type, DrugConc, sep = '_') %>% 
  select(ID, MetID, Mean) %>%
  pivot_wider(names_from = ID, values_from = Mean)

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

# annotation layer
ha = columnAnnotation(strain = c(rep('rcdA',4),rep('OP50',4)),
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




