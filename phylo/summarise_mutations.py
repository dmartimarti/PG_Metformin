#!/usr/bin/env python3

'''
This scritp goes into every folder and subfolder in the
results from snippy and merges all the info into a csv
with information about the strain used and the selected
colony from the experiment. 

Example of use, 
go to path  D://MRC_Postdoc//Pangenomic//mutants_analysis//results>

python summarise_mutations.py
'''

import os
import pandas as pd

# initialise variables 
res = pd.DataFrame()
path = os.getcwd()
folders = os.listdir(path)

out_file = 'mutations_summary.csv'

# just to know 
fold_n = 0
subfold_n = 0

# main loop
for folder in folders:
    if os.path.isdir(folder) and ".ipyn" not in folder:
        subfolders = os.listdir(folder)
        fold_n += 1
        for subfolder in subfolders:
            subfold_n += 1
            try:
                temp_df = pd.read_csv(folder+'/'+subfolder+'/snps.csv')
            except FileNotFoundError:
                continue
            temp_df.insert(0, 'Strain', folder)
            temp_df.insert(1, 'Colony', subfolder)
            res = res.append(temp_df)

print(f'\nI found {fold_n} strains, with {subfold_n} colonies in total!\n')
print(f'Saving summary into file {out_file}')
res.to_csv(out_file, index=False)