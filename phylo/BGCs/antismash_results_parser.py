#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 16:50:41 2020

@author: dani
"""

from Bio import SeqIO
import os
import pandas as pd
from optparse import OptionParser 
import re


# create a OptionParser 
# class object 
parser = OptionParser() 
  
parser.add_option("-o", "--out", 
                  dest = "output", 
                  help = "output file",  
                  metavar = "OUTPUT") 

(options, args) = parser.parse_args()

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

# retrieve txt files
def txt_parser(folder):
    path = os.getcwd()
    path_to_folder = os.path.join(path, folder,'knownclusterblast')
    files = []
    # r=root, d=directories, f = files
    for r, d, f in os.walk(path_to_folder):
        for file in f:
            if '.txt' in file :
                files.append(os.path.join(r, file))
    return(files)

def read_gbk(gbk):
    recs = [rec for rec in SeqIO.parse(gbk, "genbank")]
    return(recs[0])


def feats_extract(gbk):
    feats = [feat for feat in gbk.features if feat.type == "cand_cluster"]
    prod = feats[0].qualifiers['product']
    kind = feats[0].qualifiers['kind']
    cont = feats[0].qualifiers['contig_edge']
    if 'SMILES' in feats[0].qualifiers:
        smiles = feats[0].qualifiers['SMILES']
    else:
        smiles = 'None'
    my_dict = {'product':prod, 'kind':kind, 'contig_edge':cont, 'smiles':smiles}
    return(my_dict)


# txt parsing functions
rx_dict = {
    'source': re.compile(r'Source: (?P<source>.*)\n'),
    'type': re.compile(r'Type: (?P<type>.*)\n'),
    'blast_score': re.compile(r'Cumulative BLAST score: (?P<blast_score>\d+)\n')
}


def _parse_line(line):
    """
    Do a regex search against all defined regexes and
    return the key and match result of the first matching regex

    """

    for key, rx in rx_dict.items():
        match = rx.search(line)
        if match:
            return key, match
    # if there are no matches
    return None, None



def parse_file(filepath):
    """
    Parse text at given filepath

    Parameters
    ----------
    filepath : str
        Filepath for file_object to be parsed

    Returns
    -------
    data : pd.DataFrame
        Parsed data

    """

    data = []  # create an empty list to collect the data
    # open the file and read through it line by line
    with open(filepath, 'r') as file_object:
        line = file_object.readline()
        while line:
            # at each line check for a match with a regex
            key, match = _parse_line(line)

            # extract source
            if key == 'source':
                source = match.group('source')
            
            # extract type
            if key == 'type':
                _type = match.group('type')

            # extract blast score
            if key == 'blast_score':
                blast_score = match.group('blast_score')
                blast_score = int(blast_score)

                while line.strip():
                    # create a dictionary containing this row of data
                    row = {
                        'Source': source,
                        'Type': _type ,
                        'Blast score': blast_score
                    }
                    # append the dictionary to the data list
                    data.append(row)
                    line = file_object.readline()

            line = file_object.readline()
            newdict={}
            for k,v in [(key,d[key]) for d in data for key in d]:
                if k not in newdict: newdict[k]=[v]
                else: newdict[k].append(v)


    data = pd.DataFrame.from_dict(newdict)

        # create a pandas DataFrame from the list of dicts
        # data = pd.DataFrame(data)
        # set the School, Grade, and Student number as the index
    return data



genomes = folder_list()

# all info will be stored as diff keys in this dict

global_dict = {}
global_dict_ext = pd.DataFrame()

# read and modify gbk files
for genome in genomes:
    print('parsing genome ' + str(genome))
    list_of_files = file_parser(genome)
    list_of_txt = txt_parser(genome)
    for i in range(len(list_of_files)):
        # parse file gbk file name
        gbk_name = list_of_files[i].split('/')[-1].split('.')[0].split('_')[0]
        temp = read_gbk(list_of_files[i])
        temp_dict = feats_extract(temp)
        temp_dict['Genome'] = str(genome)
        temp_dict['Region'] = str(i+1)
        
        # update the big one
        for j in range(len(list_of_txt)):
            # parse txt file name
            clust_name = list_of_txt[j].split('/')[-1].split('.')[0].split('_')[0]
            if clust_name == gbk_name:
                ext_info = parse_file(list_of_txt[j])
                ext_info['Genome'] = genome
                ext_info['Region'] = str(i+1)
                ext_info['product'] = temp_dict['product'][0]
                ext_info['kind'] = temp_dict['kind'][0]
                ext_info['contig_edge'] = temp_dict['contig_edge'][0]
                global_dict_ext = global_dict_ext.append(ext_info)
        
        #update the global dict
        global_dict[str(genome)+ '_' +str(i)] = temp_dict

        

final_dataframe = pd.DataFrame.from_dict(global_dict, orient='index')




final_dataframe.to_csv(path + '/' + str(options.output))

global_dict_ext.to_csv(path + '/Extended_' + str(options.output))
