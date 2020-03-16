# script to analyse the metformin resistance screening with KEIO mutants

# libraries
library(tidyverse)
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(FactoMineR) # for PCA
library(factoextra) # for PCA
library(here)

options(width = 220)

# "/Users/dmarti14/Documents/MRC_Postdoc/Projects/Simren/Metformin_resistance"


# load the three replicates
rep1 = read_xlsx('HTS/191025_HTS_Keio_complete.xlsx', sheet = 'replicate1') %>%
    mutate(Drug = 'Metformin') %>%
    filter(!gene %in% c('0', 'WT', 'present'))

rep2 = read_xlsx('HTS/191025_HTS_Keio_complete.xlsx', sheet = 'replicate2') %>%
    mutate(Drug = 'Metformin') %>%
    filter(!gene %in% c('0', 'WT', 'present'))

rep3 = read_xlsx('HTS/191025_HTS_Keio_complete.xlsx', sheet = 'replicate3') %>%
    mutate(Drug = 'Metformin') %>%
    filter(!gene %in% c('0', 'WT', 'present'))

keio = rbind(rep1, rep2, rep3)

# make metformin, replicate and time factors for R interpretation
keio = keio %>%
	mutate_at(c('metformin', 'replicate', 'time'), as.factor)

## checking the data

# are the replicates the same?
# if they are, we should always have the same numbers in N (for number of cases)
keio %>% group_by(replicate) %>% summarise(N = n())
keio %>% group_by(replicate, time) %>% summarise(N = n())
keio %>% group_by(replicate, time, metformin) %>% summarise(N = n())




# check for duplicates per replicate and condition
# extract gene names
test0 = keio %>%
	filter(replicate == 1,
		   metformin == 0,
		   time == 8) %>%
	select(gene) %>% t %>% as.vector

test100 = keio %>%
	filter(replicate == 1,
		   metformin == 100,
		   time == 8) %>%
	select(gene) %>% t %>% as.vector

# these genes are duplicated
test0[duplicated(test0)]
test100[duplicated(test100)]

# are the genes the same between conditions
test0[duplicated(test0)] == test100[duplicated(test100)]

dup.genes = test0[duplicated(test0)]


# let's plot per replicate and time 



# keio %>% 
# 	filter(replicate == 1, !gene %in% dup.genes, time == 8) %>%
# 	select(gene, metformin, Drug, value) %>%
#     unite(Supp, Drug, metformin) %>%
#     spread(Supp, value) %>%
# 	ggplot(aes(x = Metformin_0, y = Metformin_100)) +
#     geom_smooth(method = 'lm') +
#     geom_abline(intercept = 0, slope = 1, colour = 'grey50') + 
#         geom_point(size = 3, alpha = .4) +
#         # geom_text_repel(aes(label = ifelse(Metformin_100 > 0.7, as.character(gene),''))) +
#         theme_classic()




require(ggpmisc)
rep = 3
ti = 8
formula = y ~ x
keio %>% 
	filter(replicate == rep, !gene %in% dup.genes, time == ti) %>%
	select(gene, metformin, Drug, value) %>%
    unite(Supp, Drug, metformin) %>%
    spread(Supp, value) %>%
	ggplot(aes(x = Metformin_0, y = Metformin_100)) +
    geom_smooth(method = 'lm') +
    xlim(0, 1) +
    ylim(0, 1) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., sep =  '~~~')), formula = formula, parse = TRUE) +
    geom_abline(intercept = 0, slope = 1, colour = 'grey50') + 
        geom_point(size = 3, alpha = .4) +
        # geom_text_repel(aes(label = ifelse(Metformin_100 > 0.7, as.character(gene),''))) +
        theme_classic()


quartz.save(file = here('exploration', paste0('Scatter_replicate',rep,'_',ti,'hours.pdf')),
    type = 'pdf', dpi = 300, height = 6, width = 8)



# # plot with linear fit table
# formula = y ~ x
# keio %>% 
# 	filter(replicate == 1, !gene %in% dup.genes, time == 8) %>%
# 	select(gene, metformin, Drug, value) %>%
#     unite(Supp, Drug, metformin) %>%
#     spread(Supp, value) %>%
# 	ggplot(aes(x = Metformin_0, y = Metformin_100)) +
#     geom_smooth(method = 'lm') +
#     stat_fit_tb(method = 'lm',
#     			method.args = list(formula = formula),
#     			tb.vars = c(Parameter = 'term',
#     						Estimate = 'estimate',
#     						"s.e." = "std.error",
#     						"italic(t)" = "statistic",
#     						"italic(P)" = "p.value"),
#     			label.y = 'top', label.x = 'left', parse = TRUE) +
#     geom_abline(intercept = 0, slope = 1, colour = 'grey50') + 
#         geom_point(size = 3, alpha = .4) +
#         # geom_text_repel(aes(label = ifelse(Metformin_100 > 0.7, as.character(gene),''))) +
#         theme_classic()




keio.sum = keio %>% group_by(time, gene, metformin, Drug) %>%
	summarise(Mean = mean(value, na.rm = T),
			  SD = sd(value, na.rm = T))





require(ggpmisc)
rep = 3
ti = 8
formula = y ~ x
keio %>% 
	filter(replicate == rep, !gene %in% dup.genes, time == ti) %>%
	select(gene, metformin, Drug, value) %>%
    unite(Supp, Drug, metformin) %>%
    spread(Supp, value) %>%
	ggplot(aes(x = Metformin_0, y = Metformin_100)) +
    geom_smooth(method = 'lm') +
    xlim(0, 1) +
    ylim(0, 1) +
    stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., sep =  '~~~')), formula = formula, parse = TRUE) +
    geom_abline(intercept = 0, slope = 1, colour = 'grey50') + 
        geom_point(size = 3, alpha = .4) +
        # geom_text_repel(aes(label = ifelse(Metformin_100 > 0.7, as.character(gene),''))) +
        theme_classic()





# quick and dirty analysis for some of Filips's questions

rep1 = keio %>% 
	filter(replicate == 1, !gene %in% dup.genes, time == 8) %>%
	select(gene, metformin, Drug, value) %>%
    unite(Supp, Drug, metformin) %>%
    spread(Supp, value)


rep2 = keio %>% 
	filter(replicate == 2, !gene %in% dup.genes, time == 8) %>%
	select(gene, metformin, Drug, value) %>%
    unite(Supp, Drug, metformin) %>%
    spread(Supp, value)


rep3 = keio %>% 
	filter(replicate == 3, !gene %in% dup.genes, time == 8) %>%
	select(gene, metformin, Drug, value) %>%
    unite(Supp, Drug, metformin) %>%
    spread(Supp, value)


value = 3
df = data.frame(gene = rep2$gene, x = rep2[,value], y = rep3[,value])
df %>% 
	ggplot(aes(x = Metformin_100, y = Metformin_100.1)) +
    geom_smooth(method = 'lm') +
    stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., sep =  '~~~')), formula = formula, parse = TRUE) +
    geom_abline(intercept = 0, slope = 1, colour = 'grey50') + 
        geom_point(size = 3, alpha = .4) +
        # geom_text_repel(aes(label = ifelse(Metformin_100 > 0.7, as.character(gene),''))) +
        theme_classic()


df %>% 
	filter(Metformin_100 > 0.05, Metformin_100.1 > 0.05) %>%
	filter(Metformin_100 < 0.6, Metformin_100.1 < 0.6) %>%
	ggplot(aes(x = Metformin_100, y = Metformin_100.1)) +
    geom_smooth(method = 'lm') +
    stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., sep =  '~~~')), formula = formula, parse = TRUE) +
    geom_abline(intercept = 0, slope = 1, colour = 'grey50') + 
        geom_point(size = 3, alpha = .4) +
        # geom_text_repel(aes(label = ifelse(Metformin_100 > 0.7, as.character(gene),''))) +
        theme_classic()



value = 2
df = data.frame(gene = rep2$gene, x = rep2[,value], y = rep3[,value])
df %>% 
	ggplot(aes(x = Metformin_0, y = Metformin_0.1)) +
    geom_smooth(method = 'lm') +
    stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., sep =  '~~~')), formula = formula, parse = TRUE) +
    geom_abline(intercept = 0, slope = 1, colour = 'grey50') + 
        geom_point(size = 3, alpha = .4) +
        # geom_text_repel(aes(label = ifelse(Metformin_100 > 0.7, as.character(gene),''))) +
        theme_classic()


df %>% 
	filter(Metformin_0 > 0.05, Metformin_0.1 > 0.05) %>%
	# filter(Metformin_100 < 0.6, Metformin_100.1 < 0.6) %>%
	ggplot(aes(x = Metformin_0, y = Metformin_0.1)) +
    geom_smooth(method = 'lm') +
    stat_poly_eq(aes(label = paste(..eq.label.., ..adj.rr.label.., sep =  '~~~')), formula = formula, parse = TRUE) +
    geom_abline(intercept = 0, slope = 1, colour = 'grey50') + 
        geom_point(size = 3, alpha = .4) +
        # geom_text_repel(aes(label = ifelse(Metformin_100 > 0.7, as.character(gene),''))) +
        theme_classic()



df2 = df %>% filter(Metformin_0 > 0.05, Metformin_0.1 > 0.05)
fit = lm(Metformin_0 ~ Metformin_0.1, data = df2)






