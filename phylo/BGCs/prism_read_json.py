#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 2021
@author: Daniel Martinez

This script uses json files as input. These json files must be analyses
from PRISM secondary metabolism website
It parses the file and saves a final output csv table with info about
clusters, features, and predicted molecules.
"""

import json
import os
import pandas as pd
import numpy as np
from tqdm import tqdm

from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

from functions import file_parser,read_json,smiles2mol

# path where the jsons are saved
path = 'D:\\MRC_Postdoc\\Pangenomic\\pangenome_analysis\\ALL\\phylo_analysis\\PRISM4\\results'

# name for the output files
output_file = os.path.join(os.path.split(path)[0],'summary.csv')
output_smiles_file = os.path.join(os.path.split(path)[0],'summary_smiles.csv')


paths,files = file_parser(path)

#####
### first part of script: generate csv with general info
#####

print(f'\nI see {len(paths)} jsons, saving first general info into a csv...\n')

# initialise dataframe
res = pd.DataFrame()

# loop over jsons in the folder
for path,file in tqdm(zip(paths,files),total=len(paths)):
	# genome = file.split('.')[0]
	genome = file[:-11] # removes '.fasta.json'
	jfile = read_json(path)
	prism = jfile['prism_results']
	prism_clusters = prism['clusters']
	for cluster_idx, cluster in enumerate(prism_clusters):
	    start = cluster['start']
	    end = cluster['end']
	    families = cluster['family']
	    types = cluster['type']
	    cluster_family = '|'.join(families)
	    cluster_type = '|'.join(types)
	    contig = cluster['contig_name']
	    # append to results
	    row = pd.DataFrame({'genome': genome, 
	    	'cluster': cluster_idx, 
	        'type': cluster_type,
	        'family': cluster_family, 'contig': contig,
	        'start': start, 'end': end }, index=[0])
	    res = res.append(row)

temp_out = output_file.split("\\")[-1]
print(f'Saving file into {temp_out}.\n')
res.to_csv(output_file, index=False)

#####
### second part of script: generate csv with products info
#####

# initialise variables
nBits = 1024 # bits for morgan fp
radius = 3 # radious of morgan fp (3 is for ecfp6)
ecfp6_name = [f'Bit_{i}' for i in range(nBits)]
smiles_res = pd.DataFrame()
df_morgan = pd.DataFrame()


print('Starting second phase: detecting predicted products by PRISM...\n')

for path,file in tqdm(zip(paths,files),total=len(paths)):
	# genome = file.split('.')[0]
	genome = file[:-11] # removes '.fasta.json'
	jfile = read_json(path)
	prism = jfile['prism_results']
	prism_clusters = prism['clusters']
	for cluster_idx, cluster in enumerate(prism_clusters):
	    pathways = cluster['biosynthetic_pathways']
	    pred_smiles = [pathway['smiles'] for pathway in pathways]
	    for smiles in pred_smiles:
	        mol = None
	        try:
	            mol = smiles2mol(smiles)
	        except ValueError:
	            continue
	        mol_weight = Descriptors.ExactMolWt(mol)
	        row = pd.concat([pd.DataFrame({'genome': genome, 
	                                       'cluster': cluster_idx, 
	                                       'smiles': smiles,
	                                       'mol_weight':mol_weight}, index=[0])], axis=1)
	        smiles_res = smiles_res.append(row)

	# fingerprints part
	fingerprints = []
	for x  in smiles_res['smiles']:
	    mol = smiles2mol(x)
	    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=nBits)
	    fingerprints.append(fp)

	ecfp6_bits = [list(l) for l in fingerprints]
	df_morgan = pd.DataFrame(ecfp6_bits, columns=ecfp6_name)
	df_morgan.insert(0, 'genome', smiles_res.genome.values)
	df_morgan.insert(1, 'cluster', smiles_res.cluster.values)
	df_morgan.insert(2, 'smiles', smiles_res.smiles.values)
	df_morgan.insert(3, 'mol_weight', smiles_res.mol_weight.values)

temp_out_smiles = output_smiles_file.split("\\")[-1]
print(f'Saving file into {temp_out_smiles}.\n')
df_morgan.to_csv(output_smiles_file, index=False)


### TO DO LIST:
# -  make phases optional by if statements
# -  parallelise all this shit, for fun