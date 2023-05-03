#!/usr/bin/env Rscript

# Load required libraries
suppressMessages(library(tidyverse))
suppressMessages(library(cowplot))
suppressMessages(library(readxl))
suppressMessages(library(rstatix))
suppressMessages(library(glue))
suppressMessages(library(argparse))

theme_set(theme_cowplot(14))

# Set up argument parser
parser <- argparse::ArgumentParser(description = "Calculate statistics and plot significant variations")
parser$add_argument("-i", "--input", 
                    help = "Input file path")
parser$add_argument("-p", "--pyseer_pheno", 
                    help = "Original phenotype used for pyseer")
parser$add_argument("-w", "--worm_pheno", 
                    help = "csv file with the worm phenotype")
parser$add_argument("-o", "--output", 
                    help = "Output file path")
args <- parser$parse_args()

# create a folder if it does not exist
# Check if the folder exists
if (!dir.exists(args$output)) {
  # Create the folder
  dir.create(args$output)
}

# get the gene name for each case
temp_str = str_split(args$input, "\\/")[[1]]
temp_str = temp_str[length(temp_str)]
gene = str_split(temp_str, '\\.')[[1]][1]

print(glue("Reading gene file {gene} \n\n\n"))

# Read input file
cat(glue("Reading data from the worm phenotype and gene {gene}: \n\n"))
protein_groups <- read_csv(args$input, show_col_types = FALSE) %>% 
  arrange(pos) %>% 
  mutate(pos = factor(pos)) %>% 
  mutate(genomes = str_replace_all(genomes, "\\[", ""),
         genomes = str_replace_all(genomes, '\\]', ""),
         genomes = str_replace_all(genomes, "\\'", "")) %>% 
  rename(genome = genomes)

# if the dataset is empty, stop the script. 
if (dim(protein_groups)[1] == 0) {
  stop("Error: protein groups file is empty or does not exist.\n", call.=FALSE)
}

# Read worm phenotype data
worm <- read_csv(args$worm_pheno, show_col_types = FALSE)

pyseer <- read_delim(args$pyseer_pheno, show_col_types = FALSE, 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

# Join datasets
worm_prot <- protein_groups %>% 
  separate_rows(genome, sep = ', ') %>% 
  left_join(worm, by = 'genome') %>% 
  drop_na(ID) %>% 
  # filter by the phenotypes present in the pyseer analysis
  filter(genome %in% pyseer$IDs)

# Filter the data
# 1. remove instances where there is only 1 aa
cat("Filtering data \n")
remove <- worm_prot %>% 
  count(pos, aa) %>% 
  filter(n == 1) %>% 
  select(-n)
if (dim(remove)[1] != 0 ){
  worm_prot <- worm_prot %>% 
    filter(!(pos %in% remove$pos & aa %in% remove$aa)) 
} 
# 2. remove the groups that only have 1 element, comparison can't be done
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
write_csv(stats,  glue("{args$output}/{gene}_stats.csv"))

# Plot significant variations
cat(glue("Trying to plot groups with significant variation... \n\n"))
if (TRUE %in% (stats$p < 0.05)) {
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
  
  ggsave(paste0(args$output,'/',gene, '_variations_plot.pdf'), plot = p,
         height = 12, width = 15)
} else {
  cat("No significant variations found in the dataset!\n")
}
