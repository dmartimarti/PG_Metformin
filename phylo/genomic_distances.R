library(tidyverse)
library(readr)
library(readxl)

distance_matrix = read_table2("distance_matrix.txt", 
                               col_names = FALSE, skip = 1)


strain_db <- read_excel("~/Documents/MRC_postdoc/Pangenomic/phylo/original_data/strain_db.xlsx")


distance_matrix = distance_matrix %>% as.matrix()

rownames(distance_matrix) = distance_matrix[,1]
distance_matrix = distance_matrix[,-1]

colnames(distance_matrix) = rownames(distance_matrix)
class(distance_matrix) <- "numeric"

distances = as.dist(distance_matrix)



library(corrplot)
col4 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "#7FFF7F",
                           "cyan", "#007FFF", "blue", "#00007F"))
corrplot(distance_matrix, method = "circle", is.corr = F, tl.pos='n', col = col4(10))

hist(distance_matrix)

distance_matrix[distance_matrix > 0.05]


#get strains 
distant = sort(rowMeans(distance_matrix))[sort(rowMeans(distance_matrix)) > 0.05]

strains = str_split(names(distant), pattern = '_')

cosa = c()
for (i in 1:length(strains)){
  cosa = c(cosa,strains[[i]][1])
}

cosa

strain_db %>% filter(ID %in% cosa)  %>% select(phylogroup)
