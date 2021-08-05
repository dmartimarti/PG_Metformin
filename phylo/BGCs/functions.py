#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 2021
@author: Daniel Martinez

Helper functions for the BGC project
"""

import os
import pandas as pd
import re
import json
from rdkit import Chem


def file_parser(folder, pattern = '.json'):
    """
    Retrieve json file list from a specific folder
    Outputs two lists, a complete file path for each file, and the 
    name of the file (json)
    """
    path = os.getcwd()
    path_to_folder = os.path.join(path, folder)
    files = []
    jsons = []
    # r=root, d=directories, f = files
    for r, d, f in os.walk(path_to_folder):
        for file in f:
            if pattern in file:
                files.append(os.path.join(r, file))
                jsons.append(file)

    return((files,jsons))



def read_json(json_file):
    """
    Reads a json file and outputs its content
    """
    with open(json_file) as f:
        content = json.load(f)
    return(content)


def smiles2mol(smiles):
    """
    reads a smiles, sanitizes it, and returns a rdkit mol
    """
    mol = Chem.MolFromSmiles(str(smiles), sanitize=True)
    Chem.SanitizeMol(mol)
    return mol


