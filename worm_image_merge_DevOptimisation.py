# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 16:27:32 2020

@author: Dani
"""


# packages
import os
import pandas as pd
from os.path import  join

# get path
path = os.getcwd()
# get folder names
folders = os.listdir()

# read files within a directory
def file_list(folder):
    temp = []
    for file in os.listdir(folder):
        if file.endswith(".csv"):
            temp.append(file)
    return(temp)


rep1_path = join(path,folders[0])
# rep2_path = join(path,folders[1])
# rep3_path = join(path,folders[2])
# rep4_path = join(path,folders[3])

# list all the files with txt extension
rep1 = file_list(rep1_path)
# rep2 = file_list(rep2_path)
# rep3 = file_list(rep3_path)
# rep4 = file_list(rep4_path)



#ugly loops
rep1df = pd.DataFrame()
for file in rep1:
    temp = pd.read_csv(join(rep1_path,file))
    
    pg = file.split(sep = '_')[1]
    metf = file.split(sep = '_')[0]
    imaging_program = file.split(sep = '_')[2]
    well = file.split(sep = '_')[3]
    zstack = well[well.find('z')+1:len(well)]
    well = well[1:well.find('z')]
    
    temp.insert(0,'PG',pg)
    temp.insert(1,'Metf', metf)
    temp.insert(2,'Well', well)
    temp.insert(3,'zStack', zstack)
    temp.insert(4, 'Imaging_program', imaging_program)
    # print('Appending dataframe ', str(file))
    rep1df = rep1df.append(temp, sort=False)
    
    
# add replicate number, join everything, save into csv
rep1df.insert(4,'Replicate',1)

total = pd.DataFrame()
total = total.append([rep1df], sort=False)

# save data into a single csv
total.to_csv(join(path,'worm_dev_optimisation.csv'), sep = ',')
