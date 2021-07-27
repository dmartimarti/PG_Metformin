# -*- coding: utf-8 -*-
"""
Script to download files from EcoRef data base

Author: Daniel Martinez-Martinez
"""

import os

import urllib.request
import pandas as pd
import gzip

# read table
st = pd.read_excel('strain_db.xlsx', sheet_name='strain_db', converters={'Assembly': str, 'Annotation': str})

DOWNLOAD_ROOT = 'https://evocellnet.github.io/ecoref/data'
ASSEMBLY = '/genomes'
ANNOTATION = '/annotations'


# function to uncompress gz files
def gzip_decompress(file_path):
    inp = gzip.GzipFile(file_path, 'rb')
    s = inp.read()
    inp.close()
    out = open(file_path[:-3], 'wb')  # opens a file without the '.gz' extension
    out.write(s)
    out.close()


# delete useless files
def delete(folder):
    files = []
    # r=root, d=directories, f = files
    for r, d, f in os.walk(folder):
        for file in f:
            if '.gz' in file:
                files.append(os.path.join(r, file))
    for file in files:
        os.remove(file)


# function to download genome files from ecoref db
def file_fetch_assembly(folder):
    os.makedirs(folder, exist_ok=True)
    st_filt = st[st.Assembly.notnull()]  # filter rows that contain info in Assembly column
    rows = st_filt.shape[0]
    for row in range(rows):
        ident = st_filt.iloc[row, :].ID  # extract ID from table
        var = st_filt.iloc[row, :].Assembly  # extract assembly number from table
        file = ident + '_' + var + '.fasta.gz'
        file_path = os.path.join(folder, file)
        url = DOWNLOAD_ROOT + ASSEMBLY + '/' + file
        urllib.request.urlretrieve(url, file_path)
        gzip_decompress(file_path=file_path)
        delete(folder)


# function to download annotation files from ecoref db
def file_fetch_annotation(folder):
    os.makedirs(folder, exist_ok=True)
    st_filt = st[st.Annotation.notnull()]  # filter rows that contain info in Annotation column
    rows = st_filt.shape[0]
    for row in range(rows):
        ident = st_filt.iloc[row, :].ID  # extract ID from table
        var = st_filt.iloc[row, :].Annotation  # extract assembly number from table
        file = ident + '_' + var + '.gff.gz'
        file_path = os.path.join(folder, file)
        url = DOWNLOAD_ROOT + ANNOTATION + '/' + file
        urllib.request.urlretrieve(url, file_path)
        gzip_decompress(file_path=file_path)
        delete(folder)


file_fetch_assembly('assemblies')
file_fetch_annotation('annotations')



