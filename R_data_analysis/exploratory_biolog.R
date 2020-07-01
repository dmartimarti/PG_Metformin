
# libraries
library(tidyverse)
library(readxl)
library(broom)
library(openxlsx)

# helping function for multiplot with ggplot2
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


data = read_csv('Summary.csv', quote = "\"") %>%
  rename(AUC_raw = `595nm_f_AUC`) %>% # `750nm_f_logAUC` data column is what we need for logAUC values
  mutate(Strain = as.character(Strain),
         Row = str_match_all(Well,'[:digit:]{1,}'), #Get plate row names from well name. 
         Col = str_match_all(Well,'[:alpha:]{1,}'), #Get plate column names from well name
         Row = factor(Row, levels = 1:12), #Make them categorical variable with set order them
         Col = factor(Col, levels = LETTERS[1:8]),
         Strain = as.factor(Strain)) %>% #Change Type column coding
  unite(ID, Strain, Plate, Well, remove = F) %>%
  select(Strain, ID, Replicate, Metformin_mM, Replicate, Plate, Well, Row, Col, AUC_raw)

data %>% filter(Strain != 'EMPTY')

# data = read.table('data_rawAUC_PG.txt', sep = '\t', header = TRUE)
#data = as_tibble(data) %>% 
#  unite(ID, Strain, Plate, Well, remove = F) 


met0 = data %>% filter(Metformin_mM == 0) %>%
  select(ID, Replicate, AUC_raw) %>%
  pivot_wider(names_from = Replicate, values_from = AUC_raw, names_prefix = 'rep_')

weird_st = met0 %>% filter(rep_1 < 0.05 | rep_2 < 0.05 | rep_3 < 0.05) %>%
  select(ID) %>% as.vector %>% t 


met0 %>% filter(!ID %in% weird_st) %>%
  ggplot(aes(x = rep_1, y = rep_2)) +
  geom_point() + 
  theme_classic()


met0  %>%
  ggplot(aes(x = rep_1, y = rep_2)) +
  geom_point() + 
  theme_classic()


met0 %>% ggplot(aes(x = rep_1, y = rep_3)) +
  geom_point() + 
  theme_classic()

met0 %>% ggplot(aes(x = rep_2, y = rep_3)) +
  geom_point() + 
  theme_classic()

model = lm(rep_1 ~ rep_2, data = met0)
summary(model)
model = lm(rep_1 ~ rep_3, data = met0)
summary(model)
model = lm(rep_2 ~ rep_3, data = met0)
summary(model)


# this is telling us something weird is happening
# lets remove evolution strains

straindb = read_xlsx('strain_db.xlsx', sheet = 'strain_db')

exp_str = straindb %>% filter(Broadphenotype == 'Evolutionexperiment') %>% select(ID) %>% t %>% as.vector

met0_filt = data %>% filter(Metformin_mM == 0) %>%
  filter(!Strain %in% exp_str) %>%
  select(ID, Replicate, AUC_raw) %>%
  pivot_wider(names_from = Replicate, values_from = AUC_raw, names_prefix = 'rep_')

# remove weird strains with low growth
# weird_st = met0_filt %>% filter(rep_1 < 0.05 | rep_2 < 0.05 | rep_3 < 0.05) %>%
#   select(ID) %>% as.vector %>% t 

# met0_filt = met0_filt %>% filter(!ID %in% weird_st)

met0_filt  %>%
  ggplot(aes(x = rep_1, y = rep_2)) +
  geom_point() + 
  theme_classic()


met0_filt %>% ggplot(aes(x = rep_1, y = rep_3)) +
  geom_point() + 
  theme_classic()

met0_filt %>% ggplot(aes(x = rep_2, y = rep_3)) +
  geom_point() + 
  theme_classic()

model = glance(lm(rep_1 ~ rep_2, data = met0_filt))
rep1_2 = glance(summary(model))
model = lm(rep_1 ~ rep_3, data = met0_filt)
rep1_3 = glance(summary(model))
model = lm(rep_2 ~ rep_3, data = met0_filt)
rep2_3 = glance(summary(model))

stats = data.frame()
stats = rbind(rep1_2, rep1_3, rep2_3)




# calculate pairwise lm functions
rep.lm = function(df, met = 0, thr = 0){
  # filter by values
  weird_st = df %>% filter(rep_1 <= thr | rep_2 <= thr | rep_3 <= thr) %>%
    select(ID) %>% as.vector %>% t 
  df = df %>% filter(!ID %in% weird_st)
  
  rep1_2 = glance(summary(lm(rep_1 ~ rep_2, data = df)))
  rep1_3 = glance(summary(lm(rep_1 ~ rep_3, data = df)))
  rep2_3 = glance(summary(lm(rep_2 ~ rep_3, data = df)))
  stats = data.frame()
  stats = rbind(rep1_2, rep1_3, rep2_3)
  stats['Metformin_mM'] = met
  stats['Comparison'] = c('rep1_2', 'rep1_3', 'rep2_3')
  return(stats %>% select(Comparison, Metformin_mM, everything()))
}


## CALCULATE STATS FROM FILTERED DATASETS

met0_filt = data %>% filter(Metformin_mM == 0) %>%
  filter(!Strain %in% exp_str) %>%
  select(ID, Replicate, AUC_raw) %>%
  pivot_wider(names_from = Replicate, values_from = AUC_raw, names_prefix = 'rep_')
met50_filt = data %>% filter(Metformin_mM == 50) %>%
  filter(!Strain %in% exp_str) %>%
  select(ID, Replicate, AUC_raw) %>%
  pivot_wider(names_from = Replicate, values_from = AUC_raw, names_prefix = 'rep_')
met100_filt = data %>% filter(Metformin_mM == 100) %>%
  filter(!Strain %in% exp_str) %>%
  select(ID, Replicate, AUC_raw) %>%
  pivot_wider(names_from = Replicate, values_from = AUC_raw, names_prefix = 'rep_')
met200_filt = data %>% filter(Metformin_mM == 200) %>%
  filter(!Strain %in% exp_str) %>%
  select(ID, Replicate, AUC_raw) %>%
  pivot_wider(names_from = Replicate, values_from = AUC_raw, names_prefix = 'rep_')


met0_stats    =  rep.lm(met0_filt, thr = 0, met = 0)
met50_stats =   rep.lm(met50_filt,thr = 0, met = 50)
met100_stats = rep.lm(met100_filt, thr = 0,met = 100)
met200_stats = rep.lm(met200_filt,thr = 0, met = 200)

data.filt_fits = rbind(met0_stats, met50_stats, met100_stats, met200_stats)

## CALCULATE STATS WITHOUT FILTERING

met0 = data %>% filter(Metformin_mM == 0) %>%
  # filter(!Strain %in% exp_str) %>%
  select(ID, Replicate, AUC_raw) %>%
  pivot_wider(names_from = Replicate, values_from = AUC_raw, names_prefix = 'rep_')
met50 = data %>% filter(Metformin_mM == 50) %>%
  # filter(!Strain %in% exp_str) %>%
  select(ID, Replicate, AUC_raw) %>%
  pivot_wider(names_from = Replicate, values_from = AUC_raw, names_prefix = 'rep_')
met100 = data %>% filter(Metformin_mM == 100) %>%
  # filter(!Strain %in% exp_str) %>%
  select(ID, Replicate, AUC_raw) %>%
  pivot_wider(names_from = Replicate, values_from = AUC_raw, names_prefix = 'rep_')
met200 = data %>% filter(Metformin_mM == 200) %>%
  # filter(!Strain %in% exp_str) %>%
  select(ID, Replicate, AUC_raw) %>%
  pivot_wider(names_from = Replicate, values_from = AUC_raw, names_prefix = 'rep_')


met0_stats    =  rep.lm(met0, thr = 0, met = 0)
met50_stats =   rep.lm(met50,thr = 0, met = 50)
met100_stats = rep.lm(met100, thr = 0,met = 100)
met200_stats = rep.lm(met200,thr = 0, met = 200)

data_fits = rbind(met0_stats, met50_stats, met100_stats, met200_stats)

# save PCA info
list_of_datasets = list('Raw_non0' = data_fits, 'NonEvo_non0' = data.filt_fits)

write.xlsx(list_of_datasets, 'Replicate_fits.xlsx', colNames = T, rowNames = T) 


# Plot R2 for the two datasets
data.filt_fits %>%
  mutate(Metformin_mM = as.factor(Metformin_mM)) %>%
  ggplot(aes(x = Comparison, y = adj.r.squared, colour = Metformin_mM)) +
  geom_point(size = 4) +
  theme_classic()

dev.copy2pdf(device = cairo_pdf,
             file = 'R2_nonEvo_non0.pdf',
             width = 6, height = 6, useDingbats = FALSE)

# Plot R2 for the two datasets
data_fits %>%
  mutate(Metformin_mM = as.factor(Metformin_mM)) %>%
  ggplot(aes(x = Comparison, y = adj.r.squared, colour = Metformin_mM)) +
  geom_point(size = 4) +
  theme_classic()

dev.copy2pdf(device = cairo_pdf,
             file = 'R2_Raw_non0.pdf',
             width = 6, height = 6, useDingbats = FALSE)


###
# Rep vs Rep plots, super ugly

# first non-filtered data
# met0
p1 = met0 %>% ggplot(aes(x = rep_1, y = rep_2)) +
  geom_point() +
  ggtitle("Metformin 0 mM") +
  geom_smooth(method = "lm") +
  theme_classic()

p2 = met0 %>% ggplot(aes(x = rep_1, y = rep_3)) +
  geom_point() +
  ggtitle("Metformin 0 mM") +
  geom_smooth(method = "lm") +
  theme_classic()

p3 = met0 %>% ggplot(aes(x = rep_2, y = rep_3)) +
  geom_point() +
  ggtitle("Metformin 0 mM") +
  geom_smooth(method = "lm") +
  theme_classic()

# met50
p4 = met50 %>% ggplot(aes(x = rep_1, y = rep_2)) +
  geom_point() +
  ggtitle("Metformin 50 mM") +
  geom_smooth(method = "lm") +
  theme_classic()

p5 = met50 %>% ggplot(aes(x = rep_1, y = rep_3)) +
  geom_point() +
  ggtitle("Metformin 50 mM") +
  geom_smooth(method = "lm") +
  theme_classic()

p6 = met50 %>% ggplot(aes(x = rep_2, y = rep_3)) +
  geom_point() +
  ggtitle("Metformin 50 mM") +
  geom_smooth(method = "lm") +
  theme_classic()

# met100
p7 = met100 %>% ggplot(aes(x = rep_1, y = rep_2)) +
  geom_point() +
  ggtitle("Metformin 100 mM") +
  geom_smooth(method = "lm") +
  theme_classic()

p8 = met100 %>% ggplot(aes(x = rep_1, y = rep_3)) +
  geom_point() +
  ggtitle("Metformin 100 mM") +
  geom_smooth(method = "lm") +
  theme_classic()

p9 = met100 %>% ggplot(aes(x = rep_2, y = rep_3)) +
  geom_point() +
  ggtitle("Metformin 100 mM") +
  geom_smooth(method = "lm") +
  theme_classic()

# met200
p10 = met200 %>% ggplot(aes(x = rep_1, y = rep_2)) +
  geom_point() +
  ggtitle("Metformin 200 mM") +
  geom_smooth(method = "lm") +
  theme_classic()

p11 = met200 %>% ggplot(aes(x = rep_1, y = rep_3)) +
  geom_point() +
  ggtitle("Metformin 200 mM") +
  geom_smooth(method = "lm") +
  theme_classic()

p12 = met200 %>% ggplot(aes(x = rep_2, y = rep_3)) +
  geom_point() +
  ggtitle("Metformin 200 mM") +
  geom_smooth(method = "lm") +
  theme_classic()



multiplot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, cols=4)

dev.copy2pdf(device = cairo_pdf,
             file = 'RepvsRep_Raw_non0.pdf',
             width = 10, height = 8, useDingbats = FALSE)


###
# Rep vs Rep plots, super ugly

# first non-filtered data
# met0
p1 = met0_filt %>% ggplot(aes(x = rep_1, y = rep_2)) +
  geom_point() +
  ggtitle("Metformin 0 mM") +
  geom_smooth(method = "lm") +
  theme_classic()

p2 = met0_filt %>% ggplot(aes(x = rep_1, y = rep_3)) +
  geom_point() +
  ggtitle("Metformin 0 mM") +
  geom_smooth(method = "lm") +
  theme_classic()

p3 = met0_filt %>% ggplot(aes(x = rep_2, y = rep_3)) +
  geom_point() +
  ggtitle("Metformin 0 mM") +
  geom_smooth(method = "lm") +
  theme_classic()

# met50
p4 = met50_filt %>% ggplot(aes(x = rep_1, y = rep_2)) +
  geom_point() +
  ggtitle("Metformin 50 mM") +
  geom_smooth(method = "lm") +
  theme_classic()

p5 = met50_filt %>% ggplot(aes(x = rep_1, y = rep_3)) +
  geom_point() +
  ggtitle("Metformin 50 mM") +
  geom_smooth(method = "lm") +
  theme_classic()

p6 = met50_filt %>% ggplot(aes(x = rep_2, y = rep_3)) +
  geom_point() +
  ggtitle("Metformin 50 mM") +
  geom_smooth(method = "lm") +
  theme_classic()

# met100
p7 = met100_filt %>% ggplot(aes(x = rep_1, y = rep_2)) +
  geom_point() +
  ggtitle("Metformin 100 mM") +
  geom_smooth(method = "lm") +
  theme_classic()

p8 = met100_filt %>% ggplot(aes(x = rep_1, y = rep_3)) +
  geom_point() +
  ggtitle("Metformin 100 mM") +
  geom_smooth(method = "lm") +
  theme_classic()

p9 = met100_filt %>% ggplot(aes(x = rep_2, y = rep_3)) +
  geom_point() +
  ggtitle("Metformin 100 mM") +
  geom_smooth(method = "lm") +
  theme_classic()

# met200
p10 = met200_filt %>% ggplot(aes(x = rep_1, y = rep_2)) +
  geom_point() +
  ggtitle("Metformin 200 mM") +
  geom_smooth(method = "lm") +
  theme_classic()

p11 = met200_filt %>% ggplot(aes(x = rep_1, y = rep_3)) +
  geom_point() +
  ggtitle("Metformin 200 mM") +
  geom_smooth(method = "lm") +
  theme_classic()

p12 = met200_filt %>% ggplot(aes(x = rep_2, y = rep_3)) +
  geom_point() +
  ggtitle("Metformin 200 mM") +
  geom_smooth(method = "lm") +
  theme_classic()



multiplot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, cols=4)

dev.copy2pdf(device = cairo_pdf,
             file = 'RepvsRep_NonEvo_non0.pdf',
             width = 10, height = 8, useDingbats = FALSE)




#############
### outlier detection
library(OutlierDetection)
library("factoextra")
library(dbscan)
# filter evo strains, set the table in wide 
# set threshold
thr = 0.1
met0_filt = data %>% filter(Metformin_mM == 0) %>%
  filter(!Strain %in% exp_str) %>%
  select(ID, Replicate, AUC_raw) %>%
  pivot_wider(names_from = Replicate, values_from = AUC_raw, names_prefix = 'rep_') 

# get the list of weird strains
weird_st = met0_filt %>% filter(rep_1 <= thr | rep_2 <= thr | rep_3 <= thr) %>%
  select(ID) %>% as.vector %>% t %>% as.character

# remove names
met0_filt = met0_filt %>% filter(!ID %in% weird_st)
X = met0_filt[,2:4]

# detect outliers 
outs = dens(met0_filt[,2:4],k=4, C=1, cutoff = 0.97)
met0_filt[outs$`Location of Outlier`,]
outs


# another method
depthout(met0_filt[,2:4], rnames = FALSE, cutoff = 0.05, boottimes = 100)


# another method
disp(met0_filt[,2:4], cutoff = 0.95, rnames = FALSE, boottimes = 1000)

# Mahalanobis
maha(met0_filt[,2:4], cutoff = 0.95, rnames = FALSE)

# k nearest
nn(met0_filt[,2:4], k = 4, cutoff = 0.98, Method = "euclidean",
   rnames = FALSE, boottimes = 100)
# kth nearest
nnk(met0_filt[,2:4], k = 4, cutoff = 0.97, Method = "euclidean",
    rnames = FALSE, boottimes = 100)

# master method
outs = OutlierDetection(met0_filt[,2:4], k = 5, cutoff = 0.98,
                 Method = "euclidean", rnames = met0_filt[,1], dispersion = TRUE)

outs
outs$`Outlier Observations`


#dbscan
db = dbscan(met0_filt[,2:4], eps = 1,  minPts = 5)
fviz_cluster(db,met0_filt[,2:4])
fviz_cluster(db,met0_filt[,2:4], axes = c(2,3))

