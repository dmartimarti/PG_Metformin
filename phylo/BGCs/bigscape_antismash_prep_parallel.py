#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 12:55:19 2020

@author: dani
"""
from Bio import SeqIO
import os
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map  # or thread_map
import multiprocessing
import time

## initialise variables

n_processes = 16 # cores to use
path = os.getcwd() # path that contains the antiSMASH folders

# path = '/home/dani/Documents/MRC_postdoc/Pangenomic/phylo/original_data/my_annotations/antismash_results'

def folder_list(path = path):
    '''
    retrieve folder list
    '''
    folders = []
    files = os.listdir(path)
    for i in files:
        if os.path.isdir(i) and '.ipynb_checkpoints' not in i:
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


def main_func(genome):
    '''
    This function modifies the gbk files from antiSMASH
    up to version 6 and saves new files with modified
    headers, ready to use with bigscape
    '''
    list_of_files = file_parser(genome)
    # print('Parsing genome ', str(genome))
    for i in range(len(list_of_files)):
        record = read_gbk(list_of_files[i])
        record.id = genome
        record.description = 'Escherichia coli ' + genome
        record.annotations['source'] = 'Escherichia coli ' + genome
        record.annotations['organism'] = 'Escherichia coli ' + genome
        path_to_save = os.path.join(path, genome, genome)+'.region00'+str(i+1)+'.gbk'
        SeqIO.write(record, path_to_save, format = 'genbank')


def move_files(genome):
    '''
    Cut files from the antiSMASH folders and paste them in the selected
    folder for the analysis

    '''
    list_of_files = file_parser_NT(genome)
    for file in list_of_files:
        file_n = file.split('/')[-1]
        if genome in file_n and 'region' in file_n: # condition to move files
            os.rename(file, 
                      os.path.join(target_path, file_n))
        else:
            pass


def count_genomes(genome):
    global n
    if len(file_parser(genome)) > 0:
        n = n+1
    else:
        n = n
    return n

start_time = time.time()

genomes = folder_list()


if __name__ == '__main__':  
    # start pool of multiprocesses
    # pool = multiprocessing.Pool()
    # pool = multiprocessing.Pool(processes=n_processes)
    
    # read and modify gbk files
    print('\nReading a modifying gbk from antiSMASH \n')
    process_map(main_func,genomes,max_workers=n_processes) # tqdm version of map

    # move created files into the right folder to be analysed
    print('\nMoving new files to analysis folder \n')
    target_path = '/mnt/d/MRC_Postdoc/Pangenomic/pangenome_analysis/ALL/phylo_analysis/BGCs/bigscape_results/gbks'
    process_map(move_files,genomes,max_workers=n_processes)

    # calculate number of genomes with BGCs
    print(f'\nCalculating the number of genomes that have BGCs, this will take a moment... \n')
    n = 0
    for genome in genomes:
        count_genomes(genome)

    print(str(len(genomes)), 'GENOMES IN TOTAL')
    print(str(round(n/len(genomes) * 100,2))+'% genomes with BGCs!' )  

end_time = round(time.time() - start_time,3)
print(f'Process complete! This took {end_time} seconds.')
