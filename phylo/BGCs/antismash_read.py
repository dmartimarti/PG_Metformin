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
import re
from tqdm import tqdm

from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

from functions import file_parser,read_json,smiles2mol

# path where the jsons are saved
path = 'D:\\MRC_Postdoc\\Pangenomic\\pangenome_analysis\\ALL\\phylo_analysis\\BGCs\\antismash_jsons'

# name for the output files, and save it in a higher level from path
output_file = os.path.join(os.path.split(path)[0],'antismash_summary.csv')
output_smiles_file = os.path.join(os.path.split(path)[0],'antismash_summary_smiles.csv')


paths,files = file_parser(path)

#####
### first part of script: generate csv with general info
#####

print(f'\nI see {len(paths)} jsons, saving first general info into a csv!')
print('Extracting features from files\n')
# initialise dataframe
res = pd.DataFrame()

# loop over jsons in the folder
for path,file in tqdm(zip(paths,files),total=len(paths)):
	# genome = file.split('.')[0]
	genome = file[:-5] # removes '.fasta.json'
	jfile = read_json(path)
	records = jfile['records']
	for record in records:
		contig = record['id']
		for feature in record['features']:
			feature_type = feature['type']
			if feature_type == 'region':
				cluster = feature
				cluster_idx = cluster['qualifiers']['region_number']
				location = cluster['location'].split(':')
				start = re.sub('\[', '', location[0])
				end = re.sub('\]', '', location[1])
				cluster_type = '|'.join(cluster['qualifiers']['product'])
				# append to results
				row = pd.DataFrame({
					'genome': genome,'cluster': cluster_idx,
					'type': cluster_type,
					'contig': contig,
					'start': start,
					'end': end }, index=[0])
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

print('Starting second phase: detecting predicted products by antiSMASH...\n')



for path,file in tqdm(zip(paths,files),total=len(paths)):
	# genome = file.split('.')[0]
	genome = file[:-5] # removes '.fasta.json'
	jfile = read_json(path)
	records = jfile['records']
	for record in records:
		contig = record['id']
		modules = record['modules']
		for module_name in modules.keys():
			module = modules[module_name]
			if 'region_predictions' in module.keys():
				preds = module['region_predictions']
				for cluster_idx in preds.keys():
					cluster = preds[cluster_idx]
					for prediction in cluster:
						smiles = prediction['smiles']
						try:
							mol = smiles2mol(smiles)
						except ValueError:
							continue
						mol_weight = Descriptors.ExactMolWt(mol)
						row = pd.concat([pd.DataFrame({
							'genome': genome,
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