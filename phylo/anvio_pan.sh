# script for anvio-6.1

# path: /media/dani/2TB Disk/MRC_Postdoc/Pangenomic/phylo/original_data/assemblies

# Generating an anvi’o genomes storage
# as fasta files are not correctly formated, need to reformat them

anvi-script-reformat-fasta NT12001_189.fasta -o NT12001_189-fixed.fa -l 0 --simplify-names
anvi-script-reformat-fasta NT12002_188.fasta -o NT12002_188-fixed.fa -l 0 --simplify-names
anvi-script-reformat-fasta NT12003_214.fasta -o NT12003_214-fixed.fa -l 0 --simplify-names
anvi-script-reformat-fasta NT12004_22.fasta  -o NT12004_22-fixed.fa -l 0 --simplify-names
anvi-script-reformat-fasta NT12005_17.fasta  -o NT12005_17-fixed.fa -l 0 --simplify-names
anvi-script-reformat-fasta NT12006_215.fasta -o NT12006_215-fixed.fa -l 0 --simplify-names
anvi-script-reformat-fasta NT12007_60.fasta  -o NT12007_60-fixed.fa -l 0 --simplify-names
anvi-script-reformat-fasta NT12008_101.fasta -o NT12008_101-fixed.fa -l 0 --simplify-names
anvi-script-reformat-fasta NT12009_154.fasta -o NT12009_154-fixed.fa -l 0 --simplify-names
anvi-script-reformat-fasta NT12010_146.fasta -o NT12010_146-fixed.fa -l 0 --simplify-names
anvi-script-reformat-fasta NT12011_132.fasta -o NT12011_132-fixed.fa -l 0 --simplify-names
anvi-script-reformat-fasta NT12012_88.fasta  -o NT12012_88-fixed.fa -l 0 --simplify-names
anvi-script-reformat-fasta NT12013_22.fasta  -o NT12013_22-fixed.fa -l 0 --simplify-names
anvi-script-reformat-fasta NT12014_22.fasta  -o NT12014_22-fixed.fa -l 0 --simplify-names
anvi-script-reformat-fasta NT12015_11.fasta  -o NT12015_11-fixed.fa -l 0 --simplify-names
anvi-script-reformat-fasta NT12016_152.fasta -o NT12016_152-fixed.fa -l 0 --simplify-names
anvi-script-reformat-fasta NT12017_35.fasta  -o NT12017_35-fixed.fa -l 0 --simplify-names
anvi-script-reformat-fasta NT12018_216.fasta -o NT12018_216-fixed.fa -l 0 --simplify-names
anvi-script-reformat-fasta NT12019_36.fasta  -o NT12019_36-fixed.fa -l 0 --simplify-names

# remove useless files
mv NT12001_189-fixed.fa NT12001_189.fasta
mv NT12002_188-fixed.fa NT12002_188.fasta
mv NT12003_214-fixed.fa NT12003_214.fasta
mv NT12004_22-fixed.fa NT12004_22.fasta
mv NT12005_17-fixed.fa NT12005_17.fasta
mv NT12006_215-fixed.fa NT12006_215.fasta
mv NT12007_60-fixed.fa NT12007_60.fasta
mv NT12008_101-fixed.fa NT12008_101.fasta
mv NT12009_154-fixed.fa NT12009_154.fasta
mv NT12010_146-fixed.fa NT12010_146.fasta
mv NT12011_132-fixed.fa NT12011_132.fasta
mv NT12012_88-fixed.fa NT12012_88.fasta
mv NT12013_22-fixed.fa NT12013_22.fasta
mv NT12014_22-fixed.fa NT12014_22.fasta
mv NT12015_11-fixed.fa NT12015_11.fasta
mv NT12016_152-fixed.fa NT12016_152.fasta
mv NT12017_35-fixed.fa NT12017_35.fasta
mv NT12018_216-fixed.fa NT12018_216.fasta
mv NT12019_36-fixed.fa NT12019_36.fasta


# lets create the contig database
# open several terminals to deal with everything
anvi-gen-contigs-database -f NT12001_189.fasta -o NT12001_189.db -n 'PG_test_1'
anvi-gen-contigs-database -f NT12002_188.fasta -o NT12002_188.db -n 'PG_test_1'
anvi-gen-contigs-database -f NT12003_214.fasta -o NT12003_214.db -n 'PG_test_1'

anvi-gen-contigs-database -f NT12004_22.fasta -o NT12004_22.db -n 'PG_test_1'
anvi-gen-contigs-database -f NT12005_17.fasta -o NT12005_17.db -n 'PG_test_1'
anvi-gen-contigs-database -f NT12006_215.fasta -o NT12006_215.db -n 'PG_test_1'

anvi-gen-contigs-database -f NT12007_60.fasta -o NT12007_60.db -n 'PG_test_1'
anvi-gen-contigs-database -f NT12008_101.fasta -o NT12008_101.db -n 'PG_test_1'
anvi-gen-contigs-database -f NT12009_154.fasta -o NT12009_154.db -n 'PG_test_1'

anvi-gen-contigs-database -f NT12010_146.fasta -o NT12010_146.db -n 'PG_test_1'
anvi-gen-contigs-database -f NT12011_132.fasta -o NT12011_132.db -n 'PG_test_1'
anvi-gen-contigs-database -f NT12012_88.fasta -o NT12012_88.db -n 'PG_test_1'

anvi-gen-contigs-database -f NT12013_22.fasta -o NT12013_22.db -n 'PG_test_1'
anvi-gen-contigs-database -f NT12014_22.fasta -o NT12014_22.db -n 'PG_test_1'
anvi-gen-contigs-database -f NT12015_11.fasta -o NT12015_11.db -n 'PG_test_1'

anvi-gen-contigs-database -f NT12016_152.fasta -o NT12016_152.db -n 'PG_test_1'
anvi-gen-contigs-database -f NT12017_35.fasta -o NT12017_35.db -n 'PG_test_1'
anvi-gen-contigs-database -f NT12018_216.fasta -o NT12018_216.db -n 'PG_test_1'
anvi-gen-contigs-database -f NT12019_36.fasta -o NT12019_36.db -n 'PG_test_1'

# iterate over all genome databases to get the COG functions
for i in $(ls *.db)
do
	anvi-run-ncbi-cogs -c $i --num-threads 12
done

# iterate over all genome databases to do hmm analysis
for i in $(ls *.db)
do
	anvi-run-hmms -c $i --num-threads 12
done

# store every db into a single dbanvi-gen-contigs-database
anvi-gen-genomes-storage -e external_PG.txt \
                         -o PG-GENOMES.db



# compute the pangenome, might take a while (for only 19 genomes!)
anvi-pan-genome -g PG-GENOMES.db \
                --project-name "Ecoli_Pan" \
                --output-dir PG_out_test \
                --num-threads 12 \
                --minbit 0.5 \
                --mcl-inflation 10 



# display genomes
anvi-display-pan -p PG_out_test/Ecoli_Pan-PAN.db -g PG-GENOMES.db


# split genomes in the different bins you find interesting
anvi-split -p PG_out_test/Ecoli_Pan-PAN.db \
           -g PG-GENOMES.db \
           -C Core \
           -o SPLIT_PANs

# display split genome
anvi-display-pan -p SPLIT_PANs/Core/PAN.db -g PG-GENOMES.db 


# use --list-collections if you want to see available collections
anvi-get-sequences-for-gene-clusters -p PG_out_test/Ecoli_Pan-PAN.db -g PG-GENOMES.db -C Core -b Better_core --concatenate-gene-clusters -o better_core.fa