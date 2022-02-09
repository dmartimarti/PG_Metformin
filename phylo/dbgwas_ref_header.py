#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import SeqIO
import re


def finder(desc):
	"""
	Finds some elements in the header description for each fasta and
	returns a string containing the original header and the DBGWAS formated
	string to be used in that pipeline
	"""
	
	# find the general tag
	z = re.findall("NC_\d*\.\d*", desc)
	general = 'DB_GWAS_general_tag=' + z[0]
	
	# find the specific tag
	x = re.findall("gene=\w+", desc)
	x = x[0].split('=')[1]
	specific = 'DB_GWAS_specific_tag='+x
	
	# find the protein tag
	y = re.findall("protein=.*$", desc)
	y = y[0].split(']')[0]
	y = y[0:]
	protein = 'DB_GWAS_protein_tag=' + y 
	
	# mix everyting
	return f'{desc};{general};{specific};{protein}'


def main():

	records = list(SeqIO.parse("GCF_000005845.2_ASM584v2_cds_from_genomic.fna", "fasta"))

	for rec in records:
		rec.description = finder(rec.description)

	SeqIO.write(records, "EcoliK12_annotation_ref_DBGWAS_2022.fasta", "fasta")


if __name__ == '__main__':
	main()
