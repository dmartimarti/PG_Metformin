#!/usr/bin/env python3

'''
This scripts takes a single fasta file as an input and splits it into n number of fastas
saved in the output folder specified in the arguments
'''

# load biopython to read fasta files
from Bio import SeqIO
import argparse
import os

# create an argument parser
parser = argparse.ArgumentParser(description='Split a fasta file into multiple files')

# input file
parser.add_argument('-i', '--input', type=str, help='Input fasta file', required=True)
# output file
parser.add_argument('-o', '--output', type=str, help='Output folder', required=True)
# numer of fastas to create
parser.add_argument('-n', '--number', type=int, help='Number of fastas to create', required=True)

# parse arguments
args = parser.parse_args()

# create output_folder if it doesn't exist
if not os.path.exists(args.output):
    os.makedirs(args.output)


# read a fasta file and split it into multiple files
def split_fasta(input_file, output_folder, number):
    # read the fasta file
    records = list(SeqIO.parse(input_file, 'fasta'))
    # calculate the number of records per file as a list of integers that sum up to the total number of records
    records_per_file = [int(len(records) / number)] * number
    # add the remaining records to the first files
    for i in range(len(records) % number):
        records_per_file[i] += 1
    # create the output files
    for i in range(number):
        # create the output file
        output_file = os.path.join(output_folder, 'file_' + str(i+1) + '.fasta')
        # write the records to the output file
        with open(output_file, 'w') as f:
            for j in range(records_per_file[i]):
                SeqIO.write(records.pop(0), f, 'fasta')



# run the function
# print information
print(f'Input file: {args.input}')
print(f'Output folder: {args.output}')
print(f'Number of fastas to create: {args.number}')

if __name__ == '__main__':
    split_fasta(args.input, args.output, args.number)

print('Done!')