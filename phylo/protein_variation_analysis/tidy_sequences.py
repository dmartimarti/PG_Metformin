#!/usr/bin/env conda run -n panaroo python
""" 
This script reads a fasta file and will rename the sequences with the original genome name. 
If there are duplicated names, it will rename them with the original name + a number.
It also removes sequences with more than 50% of gaps 
Usage example: python rename_duplicates.py -i input.fasta -p phenotype.txt -o output.fasta -g gene_presence_absence.csv
"""
import argparse
from Bio import SeqIO
import pandas as pd
import os

argparser = argparse.ArgumentParser(
    description="Renames sequences with the original genome name and filters seqs with many gaps",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)

argparser.add_argument("-i", "--input", help="Input folder with aligned fasta files")
argparser.add_argument("-p", "--phenotype", help="Original phenotype txt file used with Pyseer")
argparser.add_argument("-o", "--output", help="Output folder for the tidied fasta files")
argparser.add_argument("-g", "--gene_pa", help="Gene presence absence file")

args = argparser.parse_args()

# function to read a fasta file and return a dictionary of the sequences
def read_fasta(fasta_file):
    """ Read a fasta file with biopython """
    try:
        records = list(SeqIO.parse(fasta_file, "fasta"))
        return records
    except:
        # print in red and bold
        print("\033[1;31mError reading fasta file\033[0m")
        exit()

def read_gene_pa(gene_pa):
    """ Read a gene_pa file with pandas from a csv file"""
    try:
        gene_pa = pd.read_csv(gene_pa, engine='pyarrow')
        return gene_pa
    except:
        print("\033[1;31mError reading gene presence absence file file\033[0m")
        exit()

# read a txt file separated by tabs with pandas
def read_txt(txt_file):
    """ Read a txt file with pandas from a csv file"""
    try:
        txt = pd.read_csv(txt_file, sep="\t")
        return txt
    except:
        print("\033[1;31mError reading txt file\033[0m")
        exit()

def get_gene_metadata(gene_pa, gene):
    """ Get gene metadata from gene_pa file given a specific gene present in the file"""
    try:
        gene_metadata = gene_pa.loc[gene_pa['Gene'] == gene]
        return gene_metadata
    except:
        print("Gene not found in gene_pa file!")
        return None

# are duplicated names?
def check_duplicates(seqs):
    """ Check if there are duplicated names """
    names = [seq.name for seq in seqs]
    if len(names) == len(set(names)):
        return False
    else:
        return True

def rename_ids(fasta, gene_metadata):
    """ Rename ids with the original genome name """
    ids = [ seq.id for seq in fasta ]
    id_dict = dict(zip(gene_metadata.iloc[0,:], gene_metadata.columns))
    for i in range(len(ids)):
        if ids[i] in id_dict:
            ids[i] = id_dict[ids[i]]
    return ids

# if there are duplicated names, rename them with the original name + a number
def rename_duplicates(seqs):
    """ Rename duplicated names with the original name + a number """
    names = [seq.name for seq in seqs]
    if len(names) == len(set(names)):
        return seqs
    else:
        new_names = []
        for name in names:
            if name not in new_names:
                new_names.append(name)
            else:
                i = 1
                while name + "_" + str(i) in new_names:
                    i += 1
                new_names.append(name + "_" + str(i))
        for i in range(len(seqs)):
            seqs[i].name = new_names[i]
            seqs[i].id = new_names[i]
        return seqs

# function to remove sequences with more than 50% of gaps
def remove_gaps(seqs):
    """ Remove sequences with more than 50% of gaps """
    new_seqs = []
    for seq in seqs:
        if seq.seq.count("-") / len(seq.seq) < 0.5:
            new_seqs.append(seq)
    print("Removed " + str(len(seqs) - len(new_seqs)) + " sequences with more than 50% of gaps")
    return new_seqs


# main function
def main():
    """ Main function """

    # create output folder if it doesn't exist
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    
    # read gene_pa and pheno files
    gene_pa = read_gene_pa(args.gene_pa)
    pheno = read_txt(args.phenotype)
    pheno_ids = pheno['IDs'].tolist()

    # list fasta files in folder args.input
    fasta_files = [f for f in os.listdir(args.input) if f.endswith(".aln")]

    for fasta_file in fasta_files:
        print("Processing " + fasta_file)
        # read input file
        seqs = read_fasta(os.path.join(args.input, fasta_file))
        gene = fasta_file.split('.')[0]

        # gene metadata
        gene_metadata = get_gene_metadata(gene_pa, gene)
        # rename ids
        ids = rename_ids(seqs, gene_metadata)

        # rename the ids in the sequences
        # print("\nRenaming sequences\n")
        for i in range(len(seqs)):
            seqs[i].id = ids[i]
            seqs[i].name = ids[i]
            seqs[i].description = ""
        
        # filter sequences that are in the pheno_ids list
        # print("Filtering sequences that are not in the phenotype file\n")
        seqs = [seq for seq in seqs if seq.id in pheno_ids]
        
        # # check if there are duplicated names
        # if check_duplicates(seqs):
        #     # rename duplicated names
        #     # print in orange color and bold
        #     print("\033[1;33mDuplicated names found, renaming...\033[0m\n")
        #     seqs = rename_duplicates(seqs)
        # else:
        #     print("\033[1;32mNo duplicated names, all good!\033[0m\n")
        #     # exit

        # remove sequences with more than 50% of gaps
        print("Removing sequences with more than 50% of gaps...")
        seqs = remove_gaps(seqs)
        out_name = f'{gene}_tidy.aln'

        # write output file in the output folder from args
        print("Writing output file " + out_name + "...")
        SeqIO.write(seqs, os.path.join(args.output, out_name), "fasta")
        # SeqIO.write(seqs, args.output, "fasta")
    
    print("\033[1;32mDone\033[0m\n")

if __name__ == "__main__":
    main()
