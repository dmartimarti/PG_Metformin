# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 13:44:46 2020

@author: Dani
"""


# packages
import os
import pandas as pd

# Release information
__version__ = '0.1'
_scriptname = 'biospa_384_format'
_verdata = 'Jun 2020'
_devflag = True

# get path
path = os.getcwd()

# list all the files with txt extension
files = []
for file in os.listdir(path):
    if file.endswith(".txt"):
        files.append(file)


for file in files:
    data = pd.read_csv(file, sep = '\t', skiprows=34, skipfooter=68,
            engine='python', index_col=0)
    data = data.dropna().transpose().round(3)
    data.to_csv('new_'+file, sep = '\t')
