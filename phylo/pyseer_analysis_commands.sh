### sh script to run pyseer
conda activate assembly
# cp input names file into panaroo folder
cd /mnt/d/MRC_Postdoc/Pangenomic/pangenome_analysis/ALL/phylo_analysis/assemblies/no_evo

mash sketch -p 12 -l input_files.txt -o mash_sketch

mash dist mash_sketch.msh mash_sketch.msh| square_mash > mash.tsv

cp mash.tsv /mnt/d/MRC_Postdoc/Pangenomic/pangenome_analysis/ALL/phylo_analysis/panaroo_results_noEVO


### pyseer analysis
# go to panaroo folder
cd /mnt/d/MRC_Postdoc/Pangenomic/pangenome_analysis/ALL/phylo_analysis/panaroo_results_noEVO

# copy files with genome names and worm phenotype
cp /mnt/d/MRC_Postdoc/Pangenomic/pangenome_analysis/ALL/worm_imaging/worm_phenotype_*.txt ./pyseer_output/original_tables/

# create phylogeny from the tree generated with IQ tree
python ./pyseer/scripts/phylogeny_distance.py --lmm core_fast_tree.treefile > pyseer_output/phylogeny_K.tsv

# chose number of dimensions
scree_plot_pyseer mash.tsv

pyseer --phenotypes ./pyseer_output/original_tables/worm_phenotype_ALL.txt --lmm \
 --pres gene_presence_absence.Rtab --similarity ./pyseer_output/phylogeny_K.tsv \
 --save-m mash_mds --max-dimensions 3 --cpu 8 --min-af 0.02 --max-af 0.98  > ./pyseer_output/results/worm_phenotype_ALL.tsv


pyseer --phenotypes ./pyseer_output/original_tables/worm_phenotype_ALL_no_biofilm.txt --lmm \
 --pres gene_presence_absence.Rtab --similarity ./pyseer_output/phylogeny_K.tsv \
 --save-m mash_mds --max-dimensions 3 --cpu 8 --min-af 0.02 --max-af 0.98  > ./pyseer_output/results/worm_phenotype_ALL_no_biofilm.tsv

 pyseer --phenotypes ./pyseer_output/original_tables/worm_phenotype_AUS.txt --lmm \
 --pres gene_presence_absence.Rtab --similarity ./pyseer_output/phylogeny_K.tsv \
 --save-m mash_mds --max-dimensions 3 --cpu 8 --min-af 0.02 --max-af 0.98  > ./pyseer_output/results/worm_phenotype_AUS.tsv

  pyseer --phenotypes ./pyseer_output/original_tables/worm_phenotype_AUS_no_biofilm.txt --lmm \
 --pres gene_presence_absence.Rtab --similarity ./pyseer_output/phylogeny_K.tsv \
 --save-m mash_mds --max-dimensions 3 --cpu 8 --min-af 0.02 --max-af 0.98  > ./pyseer_output/results/worm_phenotype_AUS_no_biofilm.tsv

 pyseer --phenotypes ./pyseer_output/original_tables/worm_phenotype_ECOREF.txt --lmm \
 --pres gene_presence_absence.Rtab --similarity ./pyseer_output/phylogeny_K.tsv \
 --save-m mash_mds --max-dimensions 3 --cpu 8 --min-af 0.02 --max-af 0.98  > ./pyseer_output/results/worm_phenotype_ECOREF.tsv

 pyseer --phenotypes ./pyseer_output/original_tables/worm_phenotype_ECOREF_no_biofilm.txt --lmm \
 --pres gene_presence_absence.Rtab --similarity ./pyseer_output/phylogeny_K.tsv \
 --save-m mash_mds --max-dimensions 3 --cpu 8 --min-af 0.02 --max-af 0.98  > ./pyseer_output/results/worm_phenotype_ECOREF_no_biofilm.tsv
