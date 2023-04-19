#!/usr/bin/env Rscript

# Load required libraries
library(tidyverse)
library(cowplot)
library(readxl)
library(rstatix)
library(glue)
library(argparse)

theme_set(theme_cowplot(14))

# Set up argument parser
parser <- argparse::ArgumentParser(description = "Calculate statistics and plot significant variations")
parser$add_argument("-i", "--input", help = "Input file path")
parser$add_argument("-p", "--pyseer_pheno", help = "Original phenotype used for pyseer")
parser$add_argument("-o", "--output", help = "Output file path")
args <- parser$parse_args()


gene = str_split(args$input, '\\.')[[1]][1]

# Read input file
cat(glue("Reading data from the worm phenotype and gene {gene}: \n\n"))
protein_groups <- read_csv(args$input, col_types = cols()) %>% 
  arrange(pos) %>% 
  mutate(pos = factor(pos)) %>% 
  mutate(genomes = str_replace_all(genomes, "\\[", ""),
         genomes = str_replace_all(genomes, '\\]', ""),
         genomes = str_replace_all(genomes, "\\'", "")) %>% 
  rename(genome = genomes)

# Read worm phenotype data
worm <- read_csv('ALL_worm_FC.csv', col_types = cols())

pyseer <- read_delim(args$pyseer_pheno, 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

# Join datasets
worm_prot <- protein_groups %>% 
  separate_rows(genome, sep = ', ') %>% 
  left_join(worm, by = 'genome') %>% 
  drop_na(ID) %>% 
  # filter by the phenotypes present in the pyseer analysis
  filter(ID %in% pyseer$IDs)

# Filter the data
cat("Filtering data \n")
remove <- worm_prot %>% 
  count(pos, aa) %>% 
  filter(n == 1) %>% 
  select(-n)
if (dim(remove)[1] != 0 ){
  worm_prot <- worm_prot %>% 
    filter(!(pos %in% remove$pos & aa %in% remove$aa)) 
} 

remove_pos <- worm_prot %>% 
  count(pos, aa) %>% 
  group_by(pos) %>% 
  count() %>% 
  filter(n == 1)

if (dim(remove_pos)[1] != 0){
  worm_prot <- worm_prot %>% 
    filter(!(pos %in% remove_pos$pos)) 
}

# Calculate stats
cat("Calculating the stats per position! \n")
stats <- worm_prot %>% 
  group_by(pos) %>% 
  t_test(Mean_FC ~ aa,
         p.adjust.method = 'none',
         var.equal = FALSE,
         detailed = TRUE) %>% 
  add_significance("p")

# Write stats to output file
write_csv(stats, args$output)

# Plot significant variations
cat(glue("Trying to plot groups with significant variation... \n\n"))
if (length(stats$p < 0.05) > 0) {
  cat("Significant variations found! Plotting them \n\n")
  sig_pos <- stats %>% filter(p < 0.05) %>% distinct(pos) %>% pull(pos)
  
  p <- worm_prot %>% 
    filter(pos %in% sig_pos) %>% 
    ggplot(aes(x = aa, y = Mean_FC, fill = aa)) +
    geom_boxplot() +
    geom_point(position = position_jitterdodge(), alpha = 0.4) +
    labs(x = 'aminoacid variation',
         y = 'Mean Fold Change in acs-2 brightness') +
    facet_wrap(~pos, scales = 'free_x')
  
  ggsave(paste0(gene, '_variations_plot.pdf'), plot = p,
         height = 12, width = 15)
} else {
  cat("No significant variations found in the dataset!\n")
}
