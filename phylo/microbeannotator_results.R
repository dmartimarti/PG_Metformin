
# libraries ---------------------------------------------------------------

library(readr)
library(tidyverse)
library(ComplexHeatmap)
library(matrixStats)
library(viridis)
library(readxl)



# read data ---------------------------------------------------------------


paths = read_delim("metabolic_summary__module_completeness.tab", 
                   delim = "\t", escape_double = FALSE, 
                   trim_ws = TRUE)



metadata =  read_excel("~/Documents/MRC_postdoc/Pangenomic/metadata/MAIN_metadata.xlsx", 
                            sheet = "metadata") %>% 
  mutate(ID = str_sub(fasta, 1, -7))


metadata 


paths_matrix = paths %>% 
  select(where(is.numeric)) %>% 
  as.matrix()

# fix colnames
fixed_colnames = colnames(paths_matrix) %>% 
  str_sub(start = 1, end = -14)


colnames(paths_matrix) = fixed_colnames
rownames(paths_matrix) = paths$name


# whole matrix ####

Heatmap(paths_matrix,
        row_names_side = "left",
        col=viridis(100),
        row_names_gp = gpar(fontsize = 6),
        name = 'Pathway completeness (%)',
        # row_names_max_width = max_text_width(
        #   rownames(reduced_matrix), 
        #   gp = gpar(fontsize = 4)
        # ),
        cluster_rows = F)


quartz.save(file = 'Big_heatmap.pdf',
            type = 'pdf', dpi = 300, height = 18, width = 20)



Heatmap(paths_matrix,
        row_names_side = "left",
        col=viridis(100),
        name = 'Pathway \ncompleteness (%)',
        show_row_names = F,
        show_column_names = F,
        # row_names_max_width = max_text_width(
        #   rownames(reduced_matrix), 
        #   gp = gpar(fontsize = 4)
        # ),
        cluster_rows = F)


quartz.save(file = 'Big_heatmap_nonames.pdf',
            type = 'pdf', dpi = 300, height = 9, width = 14)



## clean up a bit the matrix ####

# calculate std and means by row
row_std = rowSds(paths_matrix)
row_means = rowMeans(paths_matrix)

# index to filter
std_index = which(row_std != 0)
mean_index = unname(which(row_means != 0))

total_index = sort(unique(intersect(std_index, mean_index) ) )

reduced_matrix = paths_matrix[total_index,]


Heatmap(reduced_matrix,
        row_names_side = "left",
        col=viridis(100),
        row_names_gp = gpar(fontsize = 6),
        name = 'Pathway \ncompleteness (%)',
        row_names_max_width = max_text_width(
          rownames(reduced_matrix), 
          gp = gpar(fontsize = 4)
        ),
        cluster_rows = F)

quartz.save(file = 'reduced_heatmap.pdf',
            type = 'pdf', dpi = 300, height = 20, width = 15)


# heatmap with no names

Heatmap(reduced_matrix,
        row_names_side = "left",
        col=viridis(100),
        row_names_gp = gpar(fontsize = 6),
        name = 'Pathway \ncompleteness (%)',
        show_row_names = F,
        show_column_names = F,
        cluster_rows = F)

quartz.save(file = 'reduced_heatmap_nonames.pdf',
            type = 'pdf', dpi = 300, height = 9, width = 14)


## even MORE reduced ####
# index to filter

# test for std thresholds
hist(row_std)

std_index = which(row_std > 10)
mean_index = unname(which(row_means != 0))

total_index = sort(unique(intersect(std_index, mean_index) ) )

reduced_matrix = paths_matrix[total_index,]

dim(reduced_matrix)

Heatmap(reduced_matrix,
        row_names_side = "left",
        col=viridis(100),
        row_names_gp = gpar(fontsize = 10),
        name = 'Pathway \ncompleteness (%)',
        row_names_max_width = max_text_width(
          rownames(reduced_matrix), 
          gp = gpar(fontsize = 6)
        ),
        cluster_rows = F)

quartz.save(file = 'reduced_heatmap_thr10.pdf',
            type = 'pdf', dpi = 300, height = 8, width = 13)






### explore these results ####

reduced_tibble = t(reduced_matrix) %>% 
  as_tibble(rownames = 'Genomes')

# betaine production

bet_0_genomes = reduced_tibble %>% 
  filter(`Betaine biosynthesis, choline => betaine` < 100) %>% 
  pull(Genomes)
  

metadata %>% 
  drop_na(Broadphenotype) %>% 
  filter(ID %in% bet_0_genomes) %>% 
  dplyr::count(Broadphenotype) %>% 
  mutate(total = sum(n),
         prop = n/total)
  
  
metadata %>% 
  drop_na(Broadphenotype) %>% 
  # filter(ID %in% bet_0_genomes) %>% 
  dplyr::count(Broadphenotype) %>% 
  mutate(total = sum(n),
         prop = n/total)
  

bet_metadata = metadata %>% 
  drop_na(Broadphenotype) %>% 
  select(ID, Broadphenotype) %>% 
  filter(!(Broadphenotype %in% c('Unknown'))) %>% 
  filter(!(Broadphenotype %in% c('Commensal strain'))) %>% 
  mutate(betaine = case_when(ID %in% bet_0_genomes ~ 'non-producer',
                             TRUE ~ 'producer')) %>% 
  distinct(ID, .keep_all = T) 

# how 
table(bet_metadata$Broadphenotype, bet_metadata$betaine)

chisq.test(bet_metadata$Broadphenotype, bet_metadata$betaine, correct=FALSE)


# grouping paths


paths_grouped = paths %>% 
  group_by(`pathway group`) %>% 
  summarise(across(where(is.numeric), mean))

paths_grouped %>% write_csv('grouped_paths.csv')

paths_longer = paths_grouped %>% 
  pivot_longer(cols = `100.fasta.faa.ko`:`SPC_4.fasta.faa.ko`, 
               names_to = 'Genome',
               values_to = 'coverage')

paths_longer %>% filter(`pathway group` == 'Beta-Lactam biosynthesis') %>% 
  arrange(desc(coverage))


paths_matrix_gr = paths_grouped %>% 
  select(where(is.numeric)) %>% 
  as.matrix()

# fix colnames
fixed_colnames = colnames(paths_matrix_gr) %>% 
  str_sub(start = 1, end = -14)


colnames(paths_matrix_gr) = fixed_colnames
rownames(paths_matrix_gr) = paths_grouped$`pathway group`

paths_matrix_gr = paths_matrix_gr[which(rowMeans(paths_matrix_gr) > 0),]

# whole matrix ####

Heatmap(paths_matrix_gr,
        row_names_side = "left",
        col=viridis(100),
        row_names_gp = gpar(fontsize = 12),
        name = 'Pathway \ncompleteness (%)',
        show_column_names = F,
        row_names_max_width = max_text_width(
          rownames(reduced_matrix), 
          gp = gpar(fontsize = 6)
        ),
        rect_gp = gpar(col = "black", lwd = 0),
        cluster_rows = F)


quartz.save(file = 'Pathways_heatmap.pdf',
            type = 'pdf', dpi = 300, height = 9, width = 13)





