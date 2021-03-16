# this has also pyseer installed
conda activate panaroo

cd /home/dani/Documents/MRC_postdoc/Pangenomic/phylo/original_data/

# this will take a lot of time
panaroo -i ./pangn_gff/*.gff -o panaroo_results_clean --clean-mode strict -t 9 -a core 

# create tree
iqtree -s core_gene_alignment.aln -pre core_tree -nt 8 -fast -m GTR


# distances 
cd /home/dani/Documents/MRC_postdoc/Pangenomic/phylo/original_data/assemblies


mash sketch -p 12 -l input_files.txt -o mash_sketch

mash dist mash_sketch.msh mash_sketch.msh| square_mash > mash.tsv


cd ../panaroo_results_clean

cp ../assemblies/mash.tsv ./


# create phylogeny from the tree generated with IQ tree
python ./pyseer_scripts/scripts/phylogeny_distance.py --lmm core_tree.treefile > pyseer_output/phylogeny_K.tsv

pyseer --phenotypes phenotypes_filtered.txt --lmm --pres gene_presence_absence.Rtab --similarity ./pyseer_output/phylogeny_K.tsv --save-m mash_mds --max-dimensions 4 --cpu 8 --min-af 0.02 --max-af 0.98  > COGs_FC_lmm.tsv


# only with worm data
pyseer --phenotypes phenotypes_worm_filtered.txt --lmm --pres gene_presence_absence.Rtab --similarity ./pyseer_output/phylogeny_K.tsv --save-m mash_mds --max-dimensions 4 --cpu 8 --min-af 0.02 --max-af 0.98  > ./pyseer_output/COGs_worm_FC_lmm.tsv


# only with worm data
pyseer --phenotypes phenotypes_bacteria_filtered.txt --lmm --pres gene_presence_absence.Rtab --similarity ./pyseer_output/phylogeny_K.tsv --save-m mash_mds --max-dimensions 4 --cpu 8 --min-af 0.02 --max-af 0.98  > ./pyseer_output/COGs_bact_score_lmm.tsv





