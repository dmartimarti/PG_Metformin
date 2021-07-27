#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 12:55:19 2020

@author: dani
"""
from Bio import SeqIO
import os


path = os.getcwd()
# path = '/home/dani/Documents/MRC_postdoc/Pangenomic/phylo/original_data/my_annotations/antismash_results'

def folder_list(path = path):
    '''
    retrieve folder list
    '''
    folders = []
    files = os.listdir(path)
    for i in files:
        if os.path.isdir(i) and 'NT' in i:
            folders.append(i)
    #path1 = os.path.join(path,folders[0])
    return(folders)


def file_parser(folder):
    """
    retrieve file list except complete tenome
    """
    path = os.getcwd()
    path_to_folder = os.path.join(path, folder)
    files = []
    # r=root, d=directories, f = files
    for r, d, f in os.walk(path_to_folder):
        for file in f:
            if '.gbk' in file and file != folder+'.gbk' and file.split('/')[-1][:2] != 'NT':
                files.append(os.path.join(r, file))
    return(files)


# retrieve file list except complete tenome
def file_parser_NT(folder):
    '''
    This is a copy of prev function that also reads files starting with NT
    '''
    path = os.getcwd()
    path_to_folder = os.path.join(path, folder)
    files = []
    # r=root, d=directories, f = files
    for r, d, f in os.walk(path_to_folder):
        for file in f:
            if '.gbk' in file and file != folder+'.gbk':
                files.append(os.path.join(r, file))
    return(files)


def read_gbk(gbk):
    '''
    Open and read gbk files
    '''
    recs = [rec for rec in SeqIO.parse(gbk, "genbank")]
    return(recs[0])


genomes = folder_list()
      
# read and modify gbk files
for genome in genomes:
    list_of_files = file_parser(genome)
    print('Parsing genome ', str(genome))
    for i in range(len(list_of_files)):
        record = read_gbk(list_of_files[i])
        record.id = genome
        path_to_save = os.path.join(path, genome, genome)+'.region00'+str(i+1)+'.gbk'
        SeqIO.write(record, path_to_save, format = 'genbank')


# move created files into the right folder to be analysed
target_path = '/home/dani/Documents/MRC_postdoc/Pangenomic/phylo/original_data/bigscape_results/gbks'
for genome in genomes:
    list_of_files = file_parser_NT(genome)
    for file in list_of_files:
        file_n = file.split('/')[-1]
        if 'NT' in file_n and 'region' in file_n: # condition to move files
            os.rename(file, 
                      os.path.join(target_path, file_n))
        else:
            pass


n = 0
for genome in genomes:
    if len(file_parser(genome)) > 0:
        n = n+1
    else:
        n = n
print(str(len(genomes)), 'GENOMES IN TOTAL')
print(str(n/len(genomes) * 100)+'% genomes with BGCs!' )  

