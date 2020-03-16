# libraries
library(tidyverse)
library(readxl)
# library(ComplexHeatmap)
# library(circlize)
library(ggrepel)
# library(PFun)
# library(forcats)
# library(FactoMineR) # for PCA
# library(factoextra) # for PCA


# "/Users/dmarti14/Documents/MRC_Postdoc/Projects/Simren/aska_library/data"


# Get timeseries data
# It will take a while
aska = read_xlsx('191001_ASKA HTS_Dani.xlsx', sheet = 'Total') %>%
    mutate(Lib = 'ASKA',
           Drug = 'Metformin') %>%
    filter(!Gene == '0')

keio = read_xlsx('191001_Keio_HTS_Dani.xlsx', sheet = 'total') %>%
    mutate(Lib = 'KEIO',
           Drug = 'Metformin') %>%
    filter(!Gene %in% c('0', 'WT'))

# some tests
aska %>% group_by(Metformin) %>% summarise(N = n())
keio %>% group_by(Metformin) %>% summarise(N = n())


# extract the vector of genes
aska.0 = aska %>% filter(Metformin == 0) %>% select(Gene) %>% t %>% as.vector
aska.100 = aska %>% filter(Metformin == 100) %>% select(Gene) %>% t %>% as.vector

# these genes are duplicated
aska.0[duplicated(aska.0)]
aska.100[duplicated(aska.100)]

# are the same genes repeated in both conditions? (they should)
aska.0[duplicated(aska.0)] %in% aska.100[duplicated(aska.100)]


# filter table with duplicated genes
aska %>% filter(Gene %in% aska.0[duplicated(aska.0)]) %>% data.frame




# general plot
aska %>% 
    filter(!Gene %in% aska.0[duplicated(aska.0)]) %>%
    select(Gene, Metformin, Drug, Value) %>%
    unite(Supp, Drug, Metformin) %>%
    spread(Supp, Value) %>%
    ggplot(aes(x = Metformin_0, y = Metformin_100)) +
        geom_point(size = 3, alpha = .4) +
        theme_classic()

quartz.save(file = here('exploration', 'aska_scatter.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 10)



# plot with names

aska %>% 
    filter(!Gene %in% aska.0[duplicated(aska.0)]) %>%
    select(Gene, Metformin, Drug, Value) %>%
    unite(Supp, Drug, Metformin) %>%
    spread(Supp, Value) %>%
    ggplot(aes(x = Metformin_0, y = Metformin_100)) +
        geom_point(size = 3, alpha = .4) +
        theme_classic() +
        geom_text_repel(aes(label = ifelse(Metformin_100 > 0.5, as.character(Gene),'')))

quartz.save(file = here('exploration', 'aska_scatter_names.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 10)



########
## KEIO

# extract the vector of genes
keio.0 = keio %>% filter(Metformin == 0) %>% select(Gene) %>% t %>% as.vector
keio.100 = keio %>% filter(Metformin == 100) %>% select(Gene) %>% t %>% as.vector

# these genes are duplicated
keio.0[duplicated(keio.0)]
keio.100[duplicated(keio.100)]

# are the same genes repeated in both conditions? (they should)
keio.0[duplicated(keio.0)] %in% keio.100[duplicated(keio.100)]


# filter table with duplicated genes
keio %>% filter(Gene %in% keio.0[duplicated(keio.0)]) %>% data.frame




# general plot
keio %>% 
    filter(!Gene %in% keio.0[duplicated(keio.0)]) %>%
    select(Gene, Metformin, Drug, Value) %>%
    unite(Supp, Drug, Metformin) %>%
    spread(Supp, Value) %>%
    ggplot(aes(x = Metformin_0, y = Metformin_100)) +
        geom_point(size = 3, alpha = .4) +
        theme_classic()

quartz.save(file = here('exploration', 'keio_scatter.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 10)



# plot with names

keio %>% 
    filter(!Gene %in% keio.0[duplicated(keio.0)]) %>%
    select(Gene, Metformin, Drug, Value) %>%
    unite(Supp, Drug, Metformin) %>%
    spread(Supp, Value) %>%
    ggplot(aes(x = Metformin_0, y = Metformin_100)) +
        geom_point(size = 3, alpha = .4) +
        theme_classic() +
        geom_text_repel(aes(label = ifelse(Metformin_100 > 0.5, as.character(Gene),'')))

quartz.save(file = here('exploration', 'keio_scatter_names.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 10)




###############
# linear models
###############

fit = lm(Value ~ Metformin, data = aska)
summary(fit)


fit = lm(Value ~ Metformin, data = keio)





temp1 = aska %>% filter(Metformin == 0) %>%
    filter(Value > 0)

test = aska %>% filter(Metformin != 0) %>%
    rbind(temp1)

fit = lm(Value ~ Metformin, data = test)
summary(fit)




temp1 = keio %>% filter(Metformin == 0) %>%
    filter(Value > 0)

test = keio %>% filter(Metformin != 0) %>%
    rbind(temp1)

fit = lm(Value ~ Metformin, data = test)
summary(fit)



#---------------------------------#



keio = read_xlsx('191001_Keio_HTS_Dani_2.xlsx', sheet = 'total') %>%
    mutate(Lib = 'KEIO',
           Drug = 'Metformin') %>%
    filter(!Gene %in% c('0', 'WT'))




# extract the vector of genes
keio.0 = keio %>% filter(Metformin == 0) %>% select(Gene) %>% t %>% as.vector
keio.100 = keio %>% filter(Metformin == 100) %>% select(Gene) %>% t %>% as.vector

# these genes are duplicated
keio.0[duplicated(keio.0)]
keio.100[duplicated(keio.100)]

# are the same genes repeated in both conditions? (they should)
keio.0[duplicated(keio.0)] %in% keio.100[duplicated(keio.100)]


# filter table with duplicated genes
keio %>% filter(Gene %in% keio.0[duplicated(keio.0)]) %>% data.frame




# general plot
keio %>% 
    filter(!Gene %in% keio.0[duplicated(keio.0)]) %>%
    select(Gene, Metformin, Drug, Value) %>%
    unite(Supp, Drug, Metformin) %>%
    spread(Supp, Value) %>%
    ggplot(aes(x = Metformin_0, y = Metformin_100)) +
        geom_point(size = 3, alpha = .4) +
        theme_classic()

quartz.save(file = here('exploration', 'keio_scatter.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 10)



# plot with names

keio %>% 
    filter(!Gene %in% keio.0[duplicated(keio.0)]) %>%
    select(Gene, Metformin, Drug, Value) %>%
    unite(Supp, Drug, Metformin) %>%
    spread(Supp, Value) %>%
    ggplot(aes(x = Metformin_0, y = Metformin_100)) +
        geom_point(size = 3, alpha = .4) +
        theme_classic() +
        geom_text_repel(aes(label = ifelse(Metformin_100 > 0.5, as.character(Gene),'')))

quartz.save(file = here('exploration', 'keio_scatter_names.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 10)




part1 = read_xlsx('191001_Keio_HTS_Dani_2.xlsx', sheet = 'part1') %>%
    mutate(Lib = 'KEIO',
           Drug = 'Metformin',
           Batch = 'Part1') %>%
    filter(!Gene %in% c('0', 'WT'))


part2 = read_xlsx('191001_Keio_HTS_Dani_2.xlsx', sheet = 'part2') %>%
    mutate(Lib = 'KEIO',
           Drug = 'Metformin',
           Batch = 'Part2') %>%
    filter(!Gene %in% c('0', 'WT'))



keio = rbind(part1, part2)


batch = keio





# general plot
keio %>% 
    filter(!Gene %in% keio.0[duplicated(keio.0)]) %>%
    select(Gene, Metformin, Drug, Value) %>%
    unite(Supp, Drug, Metformin) %>%
    spread(Supp, Value) %>%
    left_join(keio %>% select(Gene, Batch)) %>%
    ggplot(aes(x = Metformin_0, y = Metformin_100, colour = Batch)) +
        geom_point(size = 3, alpha = .4) +
        theme_classic()


quartz.save(file = here('exploration', 'keio_scatter_Batch_effect.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 10)


fit = lm(Value ~ Metformin, data = part1)
fit

residues = resid(fit)

up = residues[residues > 0.3]
down = residues[residues < -0.5]


up = as.integer(names(up))
down = as.integer(names(down))

gns = c(up, down)

gns.1 = keio[gns,]$Gene

df = data.frame(pos = gns, Gene = gns.1, Gene2 = gns.1)


# general plot
keio %>% 
    filter(!Gene %in% keio.0[duplicated(keio.0)]) %>%
    select(Gene, Metformin, Drug, Value) %>%
    unite(Supp, Drug, Metformin) %>%
    spread(Supp, Value) %>%
    left_join(keio %>% select(Gene, Batch)) %>%
    left_join(df) %>%
    ggplot(aes(x = Metformin_0, y = Metformin_100, colour = Batch)) +
        geom_point(size = 3, alpha = .4) +
        theme_classic() +
        geom_text_repel(aes(label = Gene2), colour = 'black')








############################################
############################################
############################################
############################################
############################################
############################################


# libraries
library(tidyverse)
library(readxl)
# library(ComplexHeatmap)
# library(circlize)
library(here)
library(ggrepel)
# library(PFun)
# library(forcats)
# library(FactoMineR) # for PCA
# library(factoextra) # for PCA

# WORKING PATH: /Users/dmarti14/Documents/MRC_Postdoc/Projects/Simren/Metformin_resistance/


wd = "/Users/dmarti14/Documents/MRC_Postdoc/Projects/Simren/Metformin_resistance/run1/8_hours"

setwd(wd)



data1 = read_xlsx('191009_KEIO_8hrs_normalised.xlsx', sheet = 'Plate 1')
data2 = read_xlsx('191009_KEIO_8hrs_normalised.xlsx', sheet = 'Plate 2')
data3 = read_xlsx('191009_KEIO_8hrs_normalised.xlsx', sheet = 'Plate 3')
data4 = read_xlsx('191009_KEIO_8hrs_normalised.xlsx', sheet = 'plate 4')
data5 = read_xlsx('191009_KEIO_8hrs_normalised.xlsx', sheet = 'Sheet5')
data6 = read_xlsx('191009_KEIO_8hrs_normalised.xlsx', sheet = 'Sheet6')
data7 = read_xlsx('191009_KEIO_8hrs_normalised.xlsx', sheet = 'Sheet7')
data8 = read_xlsx('191009_KEIO_8hrs_normalised.xlsx', sheet = 'Sheet8')
data9 = read_xlsx('191009_KEIO_8hrs_normalised.xlsx', sheet = 'Sheet9')
data10 = read_xlsx('191009_KEIO_8hrs_normalised.xlsx', sheet = 'Sheet10')
data11 = read_xlsx('191009_KEIO_8hrs_normalised.xlsx', sheet = 'Sheet11')
data12 = read_xlsx('191009_KEIO_8hrs_normalised.xlsx', sheet = 'Sheet12')
data13 = read_xlsx('191009_KEIO_8hrs_normalised.xlsx', sheet = 'Sheet13')

first = rbind(data1 ,data2 ,data3 ,data4 ,data5 ,data6 ,data7 ,data8 ,data9 ,data10 ,data11 ,data12 ,data13)

first = first %>% filter(!gene %in% c('0', 'WT', 'present')) %>%
    mutate(Time = 8,
           Drug = 'Metformin')

first %>% group_by(metformin) %>% summarise(N = n())


ctr.genes = first %>% filter(metformin == 0) %>% select(gene) %>% t %>% as.vector()
met.genes = first %>% filter(metformin == 100) %>% select(gene) %>% t %>% as.vector()

setdiff(ctr.genes, met.genes)



# these genes are duplicated
ctr.genes[duplicated(ctr.genes)]
met.genes[duplicated(met.genes)]

# are the same genes repeated in both conditions? (they should)
ctr.genes[duplicated(ctr.genes)] %in% met.genes[duplicated(met.genes)]


# filter table with duplicated genes
first %>% filter(gene %in% met.genes[duplicated(met.genes)]) %>% data.frame


frs =  first %>% 
    filter(!gene %in% met.genes[duplicated(met.genes)]) %>%
    select(gene, metformin, Drug, value) %>%
    unite(Supp, Drug, metformin) %>%
    spread(Supp, value)

# general plot
first %>% 
    filter(!gene %in% met.genes[duplicated(met.genes)]) %>%
    select(gene, metformin, Drug, value) %>%
    unite(Supp, Drug, metformin) %>%
    spread(Supp, value) %>%
    ggplot(aes(x = Metformin_0, y = Metformin_100)) +
    # geom_smooth(method = 'lm', formula = value ~ metformin, data = first %>% filter(!gene %in% met.genes[duplicated(met.genes)])) +
    geom_abline(intercept = fit$coefficients[2], slope = fit$coefficients[1]) + 
    geom_abline(intercept = 0, slope = 1, colour = 'grey50') + 
        geom_point(size = 3, alpha = .4) +
        geom_text_repel(aes(label = ifelse(Metformin_100 > 0.7, as.character(gene),''))) +
        theme_classic()

quartz.save(file = here('exploration', 'aska_scatter.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 10)




### 24 hours


setwd("/Users/dmarti14/Documents/MRC_Postdoc/Projects/Simren/Metformin_resistance/run1/24_hours")






data1 = read_xlsx('191009_Keio_HTS_24HRS_normalised.xlsx', sheet = 'plate1')
data2 = read_xlsx('191009_Keio_HTS_24HRS_normalised.xlsx', sheet = 'plate2')
data3 = read_xlsx('191009_Keio_HTS_24HRS_normalised.xlsx', sheet = 'plate3')
data4 = read_xlsx('191009_Keio_HTS_24HRS_normalised.xlsx', sheet = 'plate4')
data5 = read_xlsx('191009_Keio_HTS_24HRS_normalised.xlsx', sheet = 'plate5')
data6 = read_xlsx('191009_Keio_HTS_24HRS_normalised.xlsx', sheet = 'plate6')
data7 = read_xlsx('191009_Keio_HTS_24HRS_normalised.xlsx', sheet = 'plate7')
data8 = read_xlsx('191009_Keio_HTS_24HRS_normalised.xlsx', sheet = 'plate8')
data9 = read_xlsx('191009_Keio_HTS_24HRS_normalised.xlsx', sheet = 'plate9')
data10 = read_xlsx('191009_Keio_HTS_24HRS_normalised.xlsx', sheet = 'plate10')
data11 = read_xlsx('191009_Keio_HTS_24HRS_normalised.xlsx', sheet = 'plate11')
data12 = read_xlsx('191009_Keio_HTS_24HRS_normalised.xlsx', sheet = 'plate12')
data13 = read_xlsx('191009_Keio_HTS_24HRS_normalised.xlsx', sheet = 'plate13')

second = rbind(data1 ,data2 ,data3 ,data4 ,data5 ,data6 ,data7 ,data8 ,data9 ,data10 ,data11 ,data12 ,data13)

second = second %>% filter(!gene %in% c('0', 'WT', 'present')) %>%
    mutate(Time = 24,
           Drug = 'Metformin')



second %>% group_by(metformin) %>% summarise(N = n())


ctr.genes = second %>% filter(metformin == 0) %>% select(gene) %>% t %>% as.vector()
met.genes = second %>% filter(metformin == 100) %>% select(gene) %>% t %>% as.vector()

setdiff(ctr.genes, met.genes)



# these genes are duplicated
ctr.genes[duplicated(ctr.genes)]
met.genes[duplicated(met.genes)]

# are the same genes repeated in both conditions? (they should)
ctr.genes[duplicated(ctr.genes)] %in% met.genes[duplicated(met.genes)]


# filter table with duplicated genes
second %>% filter(gene %in% met.genes[duplicated(met.genes)]) %>% data.frame

fit = lm(value ~ metformin, data = second)



sec = second %>% 
    filter(!gene %in% met.genes[duplicated(met.genes)]) %>%
    select(gene, metformin, Drug, value) %>%
    unite(Supp, Drug, metformin) %>%
    spread(Supp, value)



# general plot
second %>% 
    filter(!gene %in% met.genes[duplicated(met.genes)]) %>%
    select(gene, metformin, Drug, value) %>%
    unite(Supp, Drug, metformin) %>%
    spread(Supp, value) %>%
    ggplot(aes(x = Metformin_0, y = Metformin_100)) +
    # geom_smooth(method = 'lm', formula = value ~ metformin, data = second %>% filter(!gene %in% met.genes[duplicated(met.genes)])) +
    geom_abline(intercept = fit$coefficients[2], slope = fit$coefficients[1]) + 
    geom_abline(intercept = 0, slope = 1, colour = 'grey50') + 
        geom_point(size = 3, alpha = .4) +
        # geom_text_repel(aes(label = ifelse(Metformin_100 > 0.7, as.character(gene),''))) +
        theme_classic()

quartz.save(file = here('exploration', 'aska_scatter.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 10)







# 48 hours



setwd("/Users/dmarti14/Documents/MRC_Postdoc/Projects/Simren/Metformin_resistance/run1/48_hours")






data1 = read_xlsx('191009_Keio_HTS_48HRS_normalised.xlsx', sheet  = 'plate1')
data2 = read_xlsx('191009_Keio_HTS_48HRS_normalised.xlsx', sheet  = 'plate2')
data3 = read_xlsx('191009_Keio_HTS_48HRS_normalised.xlsx', sheet  = 'plate3')
data4 = read_xlsx('191009_Keio_HTS_48HRS_normalised.xlsx', sheet  = 'plate4')
data5 = read_xlsx('191009_Keio_HTS_48HRS_normalised.xlsx', sheet  = 'plate5')
data6 = read_xlsx('191009_Keio_HTS_48HRS_normalised.xlsx', sheet  = 'plate6')
data7 = read_xlsx('191009_Keio_HTS_48HRS_normalised.xlsx', sheet  = 'plate7')
data8 = read_xlsx('191009_Keio_HTS_48HRS_normalised.xlsx', sheet  = 'plate8')
data9 = read_xlsx('191009_Keio_HTS_48HRS_normalised.xlsx', sheet  = 'plate9')
data10 = read_xlsx('191009_Keio_HTS_48HRS_normalised.xlsx', sheet = 'plate10')
data11 = read_xlsx('191009_Keio_HTS_48HRS_normalised.xlsx', sheet = 'plate11')
data12 = read_xlsx('191009_Keio_HTS_48HRS_normalised.xlsx', sheet = 'plate12')
data13 = read_xlsx('191009_Keio_HTS_48HRS_normalised.xlsx', sheet = 'plate13')

third = rbind(data1 ,data2 ,data3 ,data4 ,data5 ,data6 ,data7 ,data8 ,data9 ,data10 ,data11 ,data12 ,data13)

third = third %>% filter(!gene %in% c('0', 'WT', 'present')) %>%
    mutate(Time = 48,
           Drug = 'Metformin')



third %>% group_by(metformin) %>% summarise(N = n())


ctr.genes = third %>% filter(metformin == 0) %>% select(gene) %>% t %>% as.vector()
met.genes = third %>% filter(metformin == 100) %>% select(gene) %>% t %>% as.vector()

setdiff(ctr.genes, met.genes)



# these genes are duplicated
ctr.genes[duplicated(ctr.genes)]
met.genes[duplicated(met.genes)]

# are the same genes repeated in both conditions? (they should)
ctr.genes[duplicated(ctr.genes)] %in% met.genes[duplicated(met.genes)]


# filter table with duplicated genes
third %>% filter(gene %in% met.genes[duplicated(met.genes)]) %>% data.frame

fit = lm(value ~ metformin, data = third)


# general plot
third %>% 
    filter(!gene %in% met.genes[duplicated(met.genes)]) %>%
    select(gene, metformin, Drug, value) %>%
    unite(Supp, Drug, metformin) %>%
    spread(Supp, value) %>%
    ggplot(aes(x = Metformin_0, y = Metformin_100)) +
    # geom_smooth(method = 'lm', formula = value ~ metformin, data = second %>% filter(!gene %in% met.genes[duplicated(met.genes)])) +
    geom_abline(intercept = fit$coefficients[2], slope = fit$coefficients[1]) + 
    geom_abline(intercept = 0, slope = 1, colour = 'grey50') + 
        geom_point(size = 3, alpha = .4) +
        # geom_text_repel(aes(label = ifelse(Metformin_100 > 0.7, as.character(gene),''))) +
        theme_classic()

quartz.save(file = here('exploration', 'aska_scatter.pdf'),
    type = 'pdf', dpi = 300, height = 8, width = 10)



# sec = second %>% 
#     filter(!gene %in% met.genes[duplicated(met.genes)]) %>%
#     select(gene, metformin, Drug, value) %>%
#     unite(Supp, Drug, metformin) %>%
#     spread(Supp, value)

# frs =  first %>% 
#     filter(!gene %in% met.genes[duplicated(met.genes)]) %>%
#     select(gene, metformin, Drug, value) %>%
#     unite(Supp, Drug, metformin) %>%
#     spread(Supp, value)



# general plot
third %>% 
    filter(!gene %in% met.genes[duplicated(met.genes)]) %>%
    select(gene, metformin, Drug, value) %>%
    unite(Supp, Drug, metformin) %>%
    spread(Supp, value) %>%
    ggplot(aes(x = Metformin_0, y = Metformin_100)) +
    # geom_smooth(method = 'lm', formula = value ~ metformin, data = second %>% filter(!gene %in% met.genes[duplicated(met.genes)])) +
    geom_abline(intercept = fit$coefficients[2], slope = fit$coefficients[1]) + 
    # geom_abline(intercept = 0, slope = 1, colour = 'grey50') + 
        geom_point(size = 3, alpha = .4) +
        geom_point(data = sec, aes(x = Metformin_0, y = Metformin_100), size = 3, alpha = .4, colour = 'blue') +
        geom_point(data = frs, aes(x = Metformin_0, y = Metformin_100), size = 3, alpha = .4, colour = 'red') +
        # geom_text_repel(aes(label = ifelse(Metformin_100 > 0.7, as.character(gene),''))) +
        theme_classic()



quartz.save(file = 'scatter_8_24_48.pdf',
    type = 'pdf', dpi = 300, height = 8, width = 10)



# join all results in a single data.frame


temp1 = frs %>% mutate(Replicate = 1,
           Time_h = 8)

temp2 = sec %>% 
    mutate(Replicate = 1,
           Time_h = 24)
temp3 = third %>% 
    filter(!gene %in% met.genes[duplicated(met.genes)]) %>%
    select(gene, metformin, Drug, value) %>%
    unite(Supp, Drug, metformin) %>%
    spread(Supp, value) %>%
    mutate(Replicate = 1,
           Time_h = 48)

## COMPLETE DATASET
keio = rbind(temp1, temp2, temp3) %>% mutate(Time_h = as.factor(Time_h))

ggplot(keio, aes(Metformin_0, Metformin_100, colour = Time_h)) +
    geom_point(size = 3, alpha = .4)



##################################
### lets load second replicate ###
##################################





wd = "/Users/dmarti14/Documents/MRC_Postdoc/Projects/Simren/Metformin_resistance/run2/8_hours"

setwd(wd)



data1 = read_xlsx('191015_Keio_HTS_8hours_normalised.xlsx', sheet  = 'plate1')
data2 = read_xlsx('191015_Keio_HTS_8hours_normalised.xlsx', sheet  = 'plate2')
data3 = read_xlsx('191015_Keio_HTS_8hours_normalised.xlsx', sheet  = 'plate3')
data4 = read_xlsx('191015_Keio_HTS_8hours_normalised.xlsx', sheet  = 'plate4')
data5 = read_xlsx('191015_Keio_HTS_8hours_normalised.xlsx', sheet  = 'plate5')
data6 = read_xlsx('191015_Keio_HTS_8hours_normalised.xlsx', sheet  = 'plate6')
data7 = read_xlsx('191015_Keio_HTS_8hours_normalised.xlsx', sheet  = 'plate7')
data8 = read_xlsx('191015_Keio_HTS_8hours_normalised.xlsx', sheet  = 'plate8')
data9 = read_xlsx('191015_Keio_HTS_8hours_normalised.xlsx', sheet  = 'plate9')
data10 = read_xlsx('191015_Keio_HTS_8hours_normalised.xlsx', sheet = 'plate10')
data11 = read_xlsx('191015_Keio_HTS_8hours_normalised.xlsx', sheet = 'plate11')
data12 = read_xlsx('191015_Keio_HTS_8hours_normalised.xlsx', sheet = 'plate12')
data13 = read_xlsx('191015_Keio_HTS_8hours_normalised.xlsx', sheet = 'plate13')

first = rbind(data1 ,data2 ,data3 ,data4 ,data5 ,data6 ,data7 ,data8 ,data9 ,data10 ,data11 ,data12 ,data13)

first = first %>% filter(!gene %in% c('0', 'WT', 'present')) %>%
    mutate(Time = 8,
           Drug = 'Metformin')

first %>% group_by(metformin) %>% summarise(N = n())


ctr.genes = first %>% filter(metformin == 0) %>% select(gene) %>% t %>% as.vector()
met.genes = first %>% filter(metformin == 100) %>% select(gene) %>% t %>% as.vector()

setdiff(ctr.genes, met.genes)



# these genes are duplicated
ctr.genes[duplicated(ctr.genes)]
met.genes[duplicated(met.genes)]

# are the same genes repeated in both conditions? (they should)
ctr.genes[duplicated(ctr.genes)] %in% met.genes[duplicated(met.genes)]


# filter table with duplicated genes
first %>% filter(gene %in% met.genes[duplicated(met.genes)]) %>% data.frame


frs =  first %>% 
    filter(!gene %in% met.genes[duplicated(met.genes)]) %>%
    select(gene, metformin, Drug, value) %>%
    unite(Supp, Drug, metformin) %>%
    spread(Supp, value)





### 24 hours


setwd("/Users/dmarti14/Documents/MRC_Postdoc/Projects/Simren/Metformin_resistance/run2/24_hours")






data1 = read_xlsx('191016_Keio_HTS_24hours_normalised.xlsx', sheet = 'plate1')
data2 = read_xlsx('191016_Keio_HTS_24hours_normalised.xlsx', sheet = 'plate2')
data3 = read_xlsx('191016_Keio_HTS_24hours_normalised.xlsx', sheet = 'plate3')
data4 = read_xlsx('191016_Keio_HTS_24hours_normalised.xlsx', sheet = 'plate4')
data5 = read_xlsx('191016_Keio_HTS_24hours_normalised.xlsx', sheet = 'plate5')
data6 = read_xlsx('191016_Keio_HTS_24hours_normalised.xlsx', sheet = 'plate6')
data7 = read_xlsx('191016_Keio_HTS_24hours_normalised.xlsx', sheet = 'plate7')
data8 = read_xlsx('191016_Keio_HTS_24hours_normalised.xlsx', sheet = 'plate8')
data9 = read_xlsx('191016_Keio_HTS_24hours_normalised.xlsx', sheet = 'plate9')
data10 = read_xlsx('191016_Keio_HTS_24hours_normalised.xlsx', sheet = 'plate10')
data11 = read_xlsx('191016_Keio_HTS_24hours_normalised.xlsx', sheet = 'plate11')
data12 = read_xlsx('191016_Keio_HTS_24hours_normalised.xlsx', sheet = 'plate12')
data13 = read_xlsx('191016_Keio_HTS_24hours_normalised.xlsx', sheet = 'plate13')

second = rbind(data1 ,data2 ,data3 ,data4 ,data5 ,data6 ,data7 ,data8 ,data9 ,data10 ,data11 ,data12 ,data13)

second = second %>% filter(!gene %in% c('0', 'WT', 'present')) %>%
    mutate(Time = 24,
           Drug = 'Metformin')



second %>% group_by(metformin) %>% summarise(N = n())


ctr.genes = second %>% filter(metformin == 0) %>% select(gene) %>% t %>% as.vector()
met.genes = second %>% filter(metformin == 100) %>% select(gene) %>% t %>% as.vector()

setdiff(ctr.genes, met.genes)



# these genes are duplicated
ctr.genes[duplicated(ctr.genes)]
met.genes[duplicated(met.genes)]

# are the same genes repeated in both conditions? (they should)
ctr.genes[duplicated(ctr.genes)] %in% met.genes[duplicated(met.genes)]


# filter table with duplicated genes
second %>% filter(gene %in% met.genes[duplicated(met.genes)]) %>% data.frame

fit = lm(value ~ metformin, data = second)



sec = second %>% 
    filter(!gene %in% met.genes[duplicated(met.genes)]) %>%
    select(gene, metformin, Drug, value) %>%
    unite(Supp, Drug, metformin) %>%
    spread(Supp, value)



# 48 hours



setwd("/Users/dmarti14/Documents/MRC_Postdoc/Projects/Simren/Metformin_resistance/run2/48_hours")






data1 = read_xlsx('191017_Keio_HTS_48hours_normalised.xlsx', sheet  = 'plate1')
data2 = read_xlsx('191017_Keio_HTS_48hours_normalised.xlsx', sheet  = 'plate2')
data3 = read_xlsx('191017_Keio_HTS_48hours_normalised.xlsx', sheet  = 'plate3')
data4 = read_xlsx('191017_Keio_HTS_48hours_normalised.xlsx', sheet  = 'plate4')
data5 = read_xlsx('191017_Keio_HTS_48hours_normalised.xlsx', sheet  = 'plate5')
data6 = read_xlsx('191017_Keio_HTS_48hours_normalised.xlsx', sheet  = 'plate6')
data7 = read_xlsx('191017_Keio_HTS_48hours_normalised.xlsx', sheet  = 'plate7')
data8 = read_xlsx('191017_Keio_HTS_48hours_normalised.xlsx', sheet  = 'plate8')
data9 = read_xlsx('191017_Keio_HTS_48hours_normalised.xlsx', sheet  = 'plate9')
data10 = read_xlsx('191017_Keio_HTS_48hours_normalised.xlsx', sheet = 'plate10')
data11 = read_xlsx('191017_Keio_HTS_48hours_normalised.xlsx', sheet = 'plate11')
data12 = read_xlsx('191017_Keio_HTS_48hours_normalised.xlsx', sheet = 'plate12')
data13 = read_xlsx('191017_Keio_HTS_48hours_normalised.xlsx', sheet = 'plate13')

third = rbind(data1 ,data2 ,data3 ,data4 ,data5 ,data6 ,data7 ,data8 ,data9 ,data10 ,data11 ,data12 ,data13)

third = third %>% filter(!gene %in% c('0', 'WT', 'present')) %>%
    mutate(Time = 48,
           Drug = 'Metformin')



third %>% group_by(metformin) %>% summarise(N = n())


ctr.genes = third %>% filter(metformin == 0) %>% select(gene) %>% t %>% as.vector()
met.genes = third %>% filter(metformin == 100) %>% select(gene) %>% t %>% as.vector()

setdiff(ctr.genes, met.genes)



# these genes are duplicated
ctr.genes[duplicated(ctr.genes)]
met.genes[duplicated(met.genes)]

# are the same genes repeated in both conditions? (they should)
ctr.genes[duplicated(ctr.genes)] %in% met.genes[duplicated(met.genes)]


# filter table with duplicated genes
third %>% filter(gene %in% met.genes[duplicated(met.genes)]) %>% data.frame

fit = lm(value ~ metformin, data = third)






# join all results in a single data.frame


temp1 = frs %>% mutate(Replicate = 2,
           Time_h = 8)

temp2 = sec %>% 
    mutate(Replicate = 2,
           Time_h = 24)
temp3 = third %>% 
    filter(!gene %in% met.genes[duplicated(met.genes)]) %>%
    select(gene, metformin, Drug, value) %>%
    unite(Supp, Drug, metformin) %>%
    spread(Supp, value) %>%
    mutate(Replicate = 2,
           Time_h = 48)

## COMPLETE DATASET
keio = rbind(keio, temp1, temp2, temp3) %>% mutate(Time_h = as.factor(Time_h))

ggplot(keio, aes(Metformin_0, Metformin_100, colour = Time_h)) +
    geom_point(size = 3, alpha = .2) +
    facet_wrap(~Replicate) +
    theme_classic()

quartz.save(file = here('exploration', 'Results_per_replicate.pdf'),
    type = 'pdf', dpi = 300, height = 6, width = 12)


# summary of results


keio.sum = keio %>% group_by(gene, Time_h) %>%
    summarise(Mean_Met_0 = mean(Metformin_0, na.rm = TRUE),
              Mean_Met_100 = mean(Metformin_100, na.rm = TRUE),
              SD_Met_0 = sd(Metformin_0, na.rm = TRUE),
              SD_Met_100 = sd(Metformin_100, na.rm = TRUE))



keio.sum %>%
    ggplot(aes(x = Mean_Met_0, y = Mean_Met_100, colour = Time_h)) +
    geom_errorbar(aes(ymin = Mean_Met_100 - SD_Met_100, ymax = Mean_Met_100 + SD_Met_100), colour = 'grey50', alpha = .3) +
    geom_errorbarh(aes(xmin = Mean_Met_0 - SD_Met_0, xmax = Mean_Met_0 + SD_Met_0), colour = 'grey50', alpha = .3) +
    geom_point(size = 2, alpha = 0.8) +
    theme_classic()

quartz.save(file = here('exploration', 'Results_summary.pdf'),
    type = 'pdf', dpi = 300, height = 6, width = 8)



time = 48
keio.sum %>%
    filter(Time_h == time) %>%
    ggplot(aes(x = Mean_Met_0, y = Mean_Met_100, colour = Time_h)) +
    geom_errorbar(aes(ymin = Mean_Met_100 - SD_Met_100, ymax = Mean_Met_100 + SD_Met_100), colour = 'grey50', alpha = .3) +
    geom_errorbarh(aes(xmin = Mean_Met_0 - SD_Met_0, xmax = Mean_Met_0 + SD_Met_0), colour = 'grey50', alpha = .3) +
    geom_point(size = 2, alpha = 0.8) +
    theme_classic()

quartz.save(file = here('exploration', paste0('Results_',time,'_hours_summary.pdf')),
    type = 'pdf', dpi = 300, height = 6, width = 8)






keio.sum %>%
    # filter(Time_h == time) %>%
    ggplot(aes(x = Mean_Met_0, y = Mean_Met_100, colour = Time_h)) +
    geom_errorbar(aes(ymin = Mean_Met_100 - SD_Met_100, ymax = Mean_Met_100 + SD_Met_100), colour = 'grey50', alpha = .3) +
    geom_errorbarh(aes(xmin = Mean_Met_0 - SD_Met_0, xmax = Mean_Met_0 + SD_Met_0), colour = 'grey50', alpha = .3) +
    geom_point(size = 2, alpha = 0.8) +
    facet_wrap(~Time_h, ncol = 2) +
    theme_classic()

quartz.save(file = here('exploration', 'Results_summary_per_time.pdf'),
    type = 'pdf', dpi = 300, height = 6, width = 8)














