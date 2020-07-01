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
library(gtools)
library(openxlsx)

bliss.sc = read_csv('summary/Bliss_scores_corrected.csv')


### function for enrichment with hypergeometric test

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
  N = length(unique(drug_class$Metabolite))
  
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


# # # # # # # # # # #
### Thresholding ####
# # # # # # # # # # #


# initiate the threshold
thrs = 5
# 
# ### BE CAREFUL AND CHOSE A SPECIFIC WAY TO SEPARATE COMPOUNDS
# classify by CI_low
bliss = bliss.sc %>%
mutate(CI_low = abs(Synergy.score) - abs(CI),
       CI_up = abs(Synergy.score) + abs(CI),
       diff = abs(Synergy.score) - abs(Most.synergistic.area.score),
       Direction = ifelse(CI_low > thrs & Synergy.score < 0, 'Antagonistic',
                          ifelse(CI_low > thrs & Synergy.score > 0, 'Synergistic', 'Neutral')),
       Drug.combination = fct_reorder(Drug.combination, desc(Synergy.score)))


# 
# # classify by Score
# bliss = bliss.sc %>%
#   mutate(CI_low = abs(Synergy.score) - abs(CI),
#          CI_up = abs(Synergy.score) + abs(CI),
#          Direction = ifelse(Synergy.score < -thrs, 'Antagonistic',
#                             ifelse(Synergy.score > thrs, 'Synergistic', 'Neutral')),
#          Drug.combination = fct_reorder(Drug.combination, desc(Synergy.score)))
# 
# 


bliss %>%
  ggplot(aes(x = Drug.combination, y = Synergy.score, colour = Direction)) +
  geom_hline(yintercept = 0, size = 1, colour = 'grey50') +
  geom_errorbar(aes(x = Drug.combination, ymin = Synergy.score - CI, ymax = Synergy.score + CI)) +
  geom_point() +
  annotate(geom = "rect", xmin = 0, xmax = Inf, ymin = -thrs, ymax = thrs, # draw rectangle
           fill = "grey50", colour = "black", alpha = 0.5) +
  theme_light() +
  scale_colour_manual(values = c('#2DB814', '#8F8F8C','#BD2924')) +
  theme(axis.text.x = element_text(angle = 45, hjust  = 1)) 

# 
# # save plot
# dev.copy2pdf(device = cairo_pdf,
#              file = here('summary', 'Bliss_ordered.pdf'),
#              width = 34, height = 10, useDingbats = FALSE)



# # # # # # # # # # 
### Enrichment ####
# # # # # # # # # #

library(readxl)
biolog = read_excel("Biolog_metabolites_UPDATED_drugClass.xlsx", 
                    sheet = "Biolog_drugs")


plates = c('PM11C', 'PM12B', 'PM13B', 'PM14A', 'PM15B', 'PM16A',
           'PM17A', 'PM18C', 'PM19', 'PM20B')

# make a table with drug classes for enrichment data
drug_target = biolog %>%
  mutate(DrugConc = str_sub(MetaboliteU, -1)) %>%
  filter(DrugConc == 1) %>%
  select(Plate, Well, Index, Metabolite, EcoCycID, KEGG_ID, Target) %>%
  separate_rows(Target, sep = ', ') %>%
  drop_na(Target) %>%
  mutate(Target = as.factor(Target)) %>%
  unite(Drug.combination, Metabolite, Plate, remove = F)

# check everything is ok
drug_target

# divide data in synergistic/antagonistic

ant = bliss %>% filter(Direction == 'Antagonistic') %>% select(Drug.combination) %>% t %>% as.vector
syn = bliss %>% filter(Direction == 'Synergistic') %>% select(Drug.combination) %>% t %>% as.vector


N = length(unique(drug_class$Metabolite))


# calculate results
mol.res = enrich(syn, ant, biolog, 'Molecule')
target.res = enrich(syn, ant, biolog, 'Target') 
process.res = enrich(syn, ant, biolog, 'Process')


sig = 0.1
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
                                          size = 7), axis.text.x = element_text(face = "bold", 
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



# save statistics list
list_of_datasets = list('Process' = process.res, 'Target' = target.res , 'Molecule' = mol.res)

write.xlsx(list_of_datasets, here('Summary', 'drug_interaction_stats.xlsx'), colNames = T, rowNames = F) 







# # # # # # # # # # # # # #
### Manual enrichment ####
# # # # # # # # # # # # # #


### calculate enrichment values
# Antagonistic 
classes = unique(drug_target$Target)
enrich.ant = c()
for (class in classes){
  class.met = drug_target %>% filter(Target == class) %>% select(Drug.combination) %>% t %>% as.vector
  m = length(class.met)
  n = N - m
  k = length(ant)
  x = length(class.met[class.met %in% ant])
  fit = phyper(q = x-1, m = m, n = n, k = k, lower.tail = FALSE)
  enrich.ant = c(enrich.ant, fit)
}
# 
# enrich.ant = p.adjust(enrich.ant, method = 'fdr')
names(enrich.ant) = classes
enrich.ant = enrich.ant[enrich.ant < 0.01]


### calculate enrichment values
# Synergistic 
classes = unique(drug_target$Target)
enrich.syn = c()
for (class in classes){
  class.met = drug_target %>% filter(Target == class) %>% select(Drug.combination) %>% t %>% as.vector
  m = length(class.met)
  n = N - m
  k = length(syn)
  x = length(class.met[class.met %in% syn])
  fit = phyper(q = x-1, m = m, n = n, k = k, lower.tail = FALSE)
  enrich.syn = c(enrich.syn, fit)
}

# enrich.syn = p.adjust(enrich.syn, method = 'fdr')
names(enrich.syn) = classes
enrich.syn = enrich.syn[enrich.syn < 0.01]

# Antagonistic
enrich.ant
# Synergistic
enrich.syn




# resp = drug_target %>% filter(Target == 'redox metabolism') %>% select(Drug.combination) %>% t %>% as.vector
# syn
# ant
# 
# sum(syn %in%  resp)
# 


# # # # # # # # # # # # # #
### Fisher exact test ####
# # # # # # # # # # # # # #


# let's try one example
class.met = drug_target %>% filter(Target == 'redox metabolism') %>% select(Drug.combination) %>% t %>% as.vector
length(class.met)
g <- 75 ## Number of submitted genes
k <- 60 ## Size of the selection, i.e. submitted genes with at least one annotation in GO biological processes
m <- length(class.met) ## Number of "marked" elements, i.e. genes associated to this biological process
N <- length(unique(drug_class$Metabolite)) ## Total number of genes with some annotation in GOTERM_BP_FAT.  
n <- N - m ## Number of "non-marked" elements, i.e. genes not associated to this biological process
x <- length(class.met[class.met %in% syn]) ## Number of "marked" elements in the selection, i.e. genes of the group of interest that are associated to this biological process

p.value <-  phyper(q=x-1, m=m, n=n, k=k, lower.tail=FALSE)



## Prepare a two-dimensional contingency table
contingency.table <- data.frame(matrix(nrow=2, ncol=2))
rownames(contingency.table) <- c("predicted.target", "non.predicted")
colnames(contingency.table) <- c("class.member", "non.member")

## Assign the values one by one to make sure we put them in the right
## place (this is not necessary, we could enter the 4 values in a
## single instruction).
contingency.table["predicted.target", "class.member"] <- x ## Number of marked genes in the selection
contingency.table["predicted.target", "non.member"] <- k - x ## Number of non-marked genes in the selection
contingency.table["non.predicted", "class.member"] <- m - x ## Number of marked genes outside of the selection
contingency.table["non.predicted", "non.member"] <- n - (k - x) ## Number of non-marked genes in the selection


print(contingency.table)

## Print marginal sums
(contingency.row.sum <- apply(contingency.table, 1, sum))

(contingency.col.sum <- apply(contingency.table, 2, sum))

## Create a contingency table with marginal sums
contingency.table.margins <- cbind(contingency.table, contingency.row.sum)
contingency.table.margins <- rbind(contingency.table.margins, apply(contingency.table.margins, 2, sum))
names(contingency.table.margins) <- c(names(contingency.table), "total")
rownames(contingency.table.margins) <- c(rownames(contingency.table), "total")
print(contingency.table.margins)


## Check the total
print(sum(contingency.table)) ## The value shoudl equal N, since every

print(N)

## Run Fisher's exact test
ftest.result <- fisher.test(x=contingency.table, alternative="greater", simulate.p.value = FALSE, B = 10000000)
print(ftest.result)

print(ftest.result$p.value) ## Print the P-value of the exact test


# # # # # # # # # # # # # # # #  
### Permutation hyper test ####
# # # # # # # # # # # # # # # #

## make up some ‘true’ data

genome = unique(bliss$Drug.combination)

myGeneList = drug_target %>% filter(Target == 'translation') %>% select(Drug.combination) %>% t %>% as.vector

myGeneSet = syn

# genome = all drugs
# myGeneSet = gene set for redox met, for example
# myGeneList = drugs that are classified within redox met

dhyperRandom <- function(myGeneList, myGeneSet, genome){
  myRandomGS <- sample(genome, size = length(myGeneSet))
  myX <- length(which(myGeneList %in% myRandomGS))
  myM <- length(myGeneList)
  myN <- length(genome) - myM
  myK <- length(myGeneList)
  return(phyper(q = myX-1, m=myM, n=myN, k=myK, lower.tail = FALSE))
}

N = 100000
pvalue = c()
for(i in 1:N){
  pvalue[i] <- dhyperRandom(myGeneList, myGeneSet, genome)
}

pval = function(values, N){
  (length(values[values > 0.1]) + 1) / (N + 1)
}

pval(pvalue, N)

# let's try one example
class.met = drug_target %>% filter(Target == 'redox metabolism') %>% select(Drug.combination) %>% t %>% as.vector
length(class.met)
g <- 75 ## Number of submitted genes
k <- 60 ## Size of the selection, i.e. submitted genes with at least one annotation in GO biological processes
m <- length(class.met) ## Number of "marked" elements, i.e. genes associated to this biological process
N <- length(unique(drug_class$Metabolite)) ## Total number of genes with some annotation in GOTERM_BP_FAT.  
n <- N - m ## Number of "non-marked" elements, i.e. genes not associated to this biological process
x <- length(class.met[class.met %in% syn]) ## Number of "marked" elements in the selection, i.e. genes of the group of interest that are associated to this biological process

p.value <-  phyper(q=x-1, m=m, n=n, k=k, lower.tail=FALSE)
