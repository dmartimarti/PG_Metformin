#!/usr/bin/env conda run -n panaroo python
""" This script will take an alignment fasta file and will analyse where the variation is in the protein.
It will also plot the variation in a barplot. 
Usage example: python variation.py input.fasta 
"""

import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import seaborn as sns
import matplotlib.pyplot as plt

# parse arguments
parser = argparse.ArgumentParser(
    description="Takes an alignment and analyse where the variation is in the protein."
    )
parser.add_argument("-i", "--input", help="Input fasta file")

args = parser.parse_args()

# function to read an alignment fasta file and return the sequences
def read_fasta(fasta_file):
    """ Read a fasta file with biopython """
    try:
        records = list(SeqIO.parse(fasta_file, "fasta"))
        return records
    except:
        # print in red and bold
        print("\033[1;31mError reading fasta file, did you input a true fasta file?\033[0m")
        # end the program
        exit()

# function to get the positions where there is variation
def get_var_positions(seqs):
    """ Get the positions where there is variation """
    positions = []
    for i in range(len(seqs[0].seq)):
        bases = [seq.seq[i] for seq in seqs]
        if len(set(bases)) > 1:
            positions.append(i)
    return positions

# get a dataframe with all the variations per position and their frequency
def get_variations(seqs, positions):
    """ Get a dataframe with all the variations per position and their frequency """
    variations = []
    for pos in positions:
        bases = [seq.seq[pos] for seq in seqs]
        for base in set(bases):
            variations.append([pos, base, bases.count(base)])
    return pd.DataFrame(variations, columns=["pos", "aa", "count"])

# filter variations within low % and high % of the sequences
def filter_variations(variations, seqs, low=0.02, high=0.98):
    """ Filter variations within low % and high % of the sequences """
    vars2keep = variations[(variations["count"] > len(seqs) * low) & (variations["count"] < len(seqs) * high)]
    var_positions = list(set(vars2keep["pos"]))
    return variations[variations['pos'].isin(var_positions)]

# add frequency column to the filtered variations dataset, Try using .loc[row_indexer,col_indexer] = value instead
def add_freq(variations, seqs):
    """ Add frequency column to the filtered variations dataset after the 'count' column """
    with pd.option_context('mode.chained_assignment', None):
        variations["freq"] = variations["count"] / len(seqs)
    return variations
 
def get_groups(variations_df_filt, seqs):
    """ From the positions in the filtered variation list, create a dataframe with groups of ids that have the same variation """
    groups = []
    for pos in set(variations_df_filt["pos"]):
        for base in set(variations_df_filt[variations_df_filt["pos"] == pos]["aa"]):
            group = []
            for seq in seqs:
                if seq.seq[pos] == base:
                    group.append(seq.id)
            groups.append([pos, base, len(group), group])
    return pd.DataFrame(groups, columns=["pos", "aa", "count", "genomes"])

def plot_variations(variations, gene):
    """ Plot the variations """
    sns.set(style="whitegrid")
    sns.set(rc={'figure.figsize':(14.7,10.27)})
    sns.set_style("whitegrid", {'axes.grid' : False})
    sns.set_context("paper", font_scale=1.5)
    sns.barplot(x="pos", y="freq", hue="aa", width=1.4, data=variations.sort_values(by="freq", ascending=False))
    # add title
    plt.title(f"Variations in the protein: {gene}")
    # add x label
    plt.xlabel("Position")
    # add y label
    plt.ylabel("Frequency")
    # save figure
    plt.savefig(f"{gene}_variations.png", dpi=300, bbox_inches='tight')

# def main and print information in each step, with colours
def main():
    """ Main function """
    # read fasta file
    print("\033[1;32mReading fasta file...\033[1;m")
    seqs = read_fasta(args.input)
    # get gene name
    gene = args.input.split("/")[-1].split(".")[0]
    # get positions where there is variation
    print("\033[1;32mGetting positions where there is variation...\033[1;m")
    positions = get_var_positions(seqs)
    # get variations dataframe
    print("\033[1;32mGetting variations dataframe...\033[1;m")
    variations = get_variations(seqs, positions)
    # filter variations
    print("\033[1;32mFiltering out variations with af below 2% and above 98%...\033[1;m")
    variations = filter_variations(variations, seqs)
    # add frequency column
    variations = add_freq(variations, seqs)
    # get groups
    print("\033[1;32mGetting groups...\033[1;m")
    groups = get_groups(variations, seqs)
    # save groups dataframe
    print("\033[1;32mSaving groups dataframe...\033[1;m")
    groups.to_csv(f"{gene}_groups.csv", index=False)
    # plot variations
    print("\033[1;32mPlotting variations...\033[1;m")
    plot_variations(variations, gene)
    # save variations dataframe
    # print("\033[1;32mSaving variations dataframe...\033[1;m")
    # variations.to_csv(f"{gene}_variations.csv", index=False)
    print("\033[1;32mDone!\033[1;m")

if __name__ == "__main__":
    main()
