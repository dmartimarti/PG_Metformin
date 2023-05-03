#!/usr/bin/env conda run -n panaroo python
"""
This script will take a txt output from pyseer and will look for the gene names in the gene_presence_absence.csv file.
Then it will run a terminal command to get all the gene sequences from the pangenome.
Usage: python get_gene_sequences.py -i input.txt -g gene_presence_absence.csv 
"""
import difflib
import multiprocessing
from tqdm import tqdm
import pandas as pd
import argparse
import os

# parse arguments
parser = argparse.ArgumentParser(
    description="Get gene sequences from pangenome",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument("-i", "--input", help="Input txt file from pyseer")
parser.add_argument("-g", "--gene_pa", help="Gene presence absence file")

args = parser.parse_args()

# function to read a txt file separated by tabs with pandas
def read_txt(txt_file):
    """ Read a txt file with pandas from a csv file"""
    try:
        txt = pd.read_csv(txt_file, sep="\t")
        return txt
    except:
        print("\033[1;31mError reading txt file\033[0m")
        exit()

# function to read a gene_pa file with pandas from a csv file
def read_gene_pa(gene_pa):
    """ Read a gene_pa file with pandas from a csv file"""
    try:
        gene_pa = pd.read_csv(gene_pa, engine='pyarrow')
        return gene_pa
    except:
        print("\033[1;31mError reading gene presence absence file file\033[0m")
        exit()

# Function to find similar genes
def find_similar_genes(args):
    query_gene, objective_genes = args
    similar_genes = []
    for objective_gene in objective_genes:
        if query_gene in objective_gene or objective_gene in query_gene:
            similar_genes.append(objective_gene)
            continue
        elif difflib.SequenceMatcher(None, query_gene, objective_gene).ratio() >= 0.8:
            similar_genes.append(objective_gene)
    return similar_genes

# main function
def main():
    # read input txt file
    pyseer_out = read_txt(args.input)
    # read gene_pa file
    gene_pa = read_gene_pa(args.gene_pa)
    # get the list of genes in the txt file
    # you can adjust the p-value threshold here!!!
    query_genes = pyseer_out[pyseer_out['maxp'] > 4].gene.unique().tolist()
    # objective gene list
    objective_genes = gene_pa['Gene']
    # Find similar genes using multiprocessing and tqdm
    print("Finding similar genes\n")
    similar_genes = []
    with multiprocessing.get_context("fork").Pool(8) as pool:
        for results in tqdm(pool.imap_unordered(find_similar_genes, [(query_gene, objective_genes) for query_gene in query_genes])):
            similar_genes.extend(results)
    # get the unique genes
    similar_genes = list(set(similar_genes))
    
    # now run the terminal command in the level above (like cd ..): 
    # panaroo-extract-gene --pa gene_presence_absence.csv --gene gene_data.csv -o genes -q gene
    # per gene
    # using multicore and tqdm
    print("Extracting genes\n")
    with multiprocessing.get_context("fork").Pool(8) as pool:
        for _ in tqdm(pool.imap_unordered(os.system, ["panaroo-extract-gene --pa ../gene_presence_absence.csv --gene ../gene_data.csv -o genes -q " + gene for gene in similar_genes])):
            pass
  
if __name__ == "__main__":
    main()