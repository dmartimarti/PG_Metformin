#!/bin/bash

if [[ ! -d antismash_dataset ]]; then
	mkdir antismash_dataset
fi


# create folder for taxonomy
if [[ ! -d taxonomy ]]; then
	mkdir taxonomy
fi

# start writting the taxonomy file
printf "# Genome folder\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tOrganism\n">> taxonomy/taxonomy_antismash.tsv

echo "Starting to get files in shape..."

for folder in $(ls ./subset/); 
# reads folders in subset and creates convenient folders for analysis
do
	echo "copying files from genome $folder \n"
	# this will print a line with dataset info into the taxonomy
	printf "$folder/\tBacteria\tProteobacteria\tGammaproteobacteria\tEnterobacteriales\tEnterobacteriaceae\tEscherichia\tEscherichia coli\tEscherichia coli $folder\n">> taxonomy/taxonomy_antismash.tsv
	# statement to create folders if they don't exist
	if [[ ! -d antismash_dataset/$folder ]]; then
		mkdir antismash_dataset/$folder # creates a folder
	fi

	# copy gbk region files into subfolders
	for file in $(ls ./subset/$folder/*region*gbk);
	do
		cp $file ./antismash_dataset/$folder
	done
done

# create dataset file with proper info
printf "# Dataset name\tPath to folder\tPath to taxonomy\tDescription\n" >> datasets.tsv
printf "antismash_dataset\tantismash_dataset/\ttaxonomy/taxonomy_antismash.tsv\tDescription" >> datasets.tsv


echo "Everything should be ready to run with BigSlice! Have fun"



