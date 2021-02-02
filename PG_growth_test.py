# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 10:30:13 2020

@author: Dani
"""


# packages
import os
import pandas as pd
from os.path import  join
import numpy as np

# get path
path = os.getcwd()
# get folder names
folders = os.listdir()

# read files within a directory
def file_list(folder):
    temp = []
    for file in os.listdir(folder):
        if file.endswith(".txt"):
            temp.append(file)
    return(temp)

# get folders paths for each time point
fold_paths = []
for folder in folders:
    fold_paths.append(join(path, folder))

# initiate variable wells
wells = ['A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12',
'B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11','B12',
'C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12',
'D1','D2','D3','D4','D5','D6','D7','D8','D9','D10','D11','D12',
'E1','E2','E3','E4','E5','E6','E7','E8','E9','E10','E11','E12',
'F1','F2','F3','F4','F5','F6','F7','F8','F9','F10','F11','F12',
'G1','G2','G3','G4','G5','G6','G7','G8','G9','G10','G11','G12',
'H1','H2','H3','H4','H5','H6','H7','H8','H9','H10','H11','H12']


full = pd.DataFrame()
for path in fold_paths:
    for file in file_list(path):
        temp = pd.read_csv(join(path,file), sep = '\t', skiprows=36, 
                   engine='python')
        # generate a clean dataset
        clean = pd.DataFrame()
        # variables that depend on the file being read
        t = path.split(sep = '\\')[-1].split(sep = '_')[0]
        plate = file[:-4].split(sep = '_')[1]
        metf = file[:-4].split(sep = '_')[2]

        # change column names
        temp.columns = [1,2,3,4,5,6,7,8,9,10,11,12,13]
        temp = temp.drop(columns = [13])
        
        # flatten the list
        flattened_list = [y for x in temp.values.tolist() for y in x]
        values = flattened_list
        
        # append variables from each file
        clean['Well'] = wells
        clean['OD']  = values 
        clean['plate'] = plate
        clean['time_h'] = int(t)
        clean['Metformin_mM'] = int(metf)
        
        full = full.append(clean, sort=False)


# arrange by values
full = full.sort_values(by = ['plate', 'time_h', 'Metformin_mM'])

# save file 
path = os.getcwd()
full.to_csv(join(path,'Summary.csv'))
