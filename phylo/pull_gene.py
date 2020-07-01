# -*- coding: utf-8 -*-
"""
Getting gene sequences from our PG files

Author: Daniel Martinez-Martinez
"""

import os

import pprint
from BCBio import GFF
from BCBio.GFF import GFFExaminer


ann_folder = '.\\annotations'
file = "NT12004_22.gff"

in_file = os.path.join(ann_folder, file)

# to explore gff files
examiner = GFFExaminer()
in_handle = open(in_file)
pprint.pprint(examiner.parent_child_map(in_handle))
# in_handle.close()

# to parse the document
in_handle = open(file)
for rec in GFF.parse(in_handle):
    print(rec)
in_handle.close()
