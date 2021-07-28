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

from functions import file_parser,read_json

path = 'D:\\MRC_Postdoc\\Pangenomic\\pangenome_analysis\\ALL\\phylo_analysis\\PRISM4\\results'

output_file = os.path.join(os.path.split(path)[0],'summary.csv')

paths,files = file_parser(path)

# initialise dataframe
res = pd.DataFrame()

# loop over 
for path,file in tqdm(zip(paths,files),total=len(paths)):
	genome = file.split('.')[0]
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

res.to_csv(output_file, index=False)