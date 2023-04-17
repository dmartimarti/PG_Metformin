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
 --save-m mash_mds --max-dimensions 3 --cpu 8 --min-af 0.02 --max-af 0.98  > ./pyseer_output/results/worm_phenotype_ALL.tsv 2>log.err


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

# Biofilm production for Jen
 pyseer --phenotypes ./pyseer_output/original_tables/worm_phenotype_normal2superbio.txt --lmm \
 --pres gene_presence_absence.Rtab --similarity ./pyseer_output/phylogeny_K.tsv \
 --save-m mash_mds --max-dimensions 3 --cpu 8 --min-af 0.02 --max-af 0.98  > ./pyseer_output/results/worm_phenotype_normal2superbio.tsv


### results with elastic nets

pyseer --phenotypes ./pyseer_output/original_tables/worm_phenotype_ALL.txt --wg enet --pres gene_presence_absence.Rtab  --lineage-clusters ./pyseer_output/phylogeny_K.tsv \
 --save-m mash_mds --max-dimensions 3 --cpu 8 --min-af 0.02 --max-af 0.98  > ./pyseer_output/results/worm_phenotype_ALL_enet.tsv 




### for bacterial growth 

# copy files with genome names and worm phenotype
cp /mnt/d/MRC_Postdoc/Pangenomic/pangenome_analysis/ALL/bacterial_growth/bact_phenotype_*.txt ./pyseer_output/original_tables/bact

pyseer --phenotypes ./pyseer_output/original_tables/bact/bact_phenotype_ALL.txt --lmm \
 --pres gene_presence_absence.Rtab --similarity ./pyseer_output/phylogeny_K.tsv \
 --save-m mash_mds --max-dimensions 3 --cpu 8 --min-af 0.02 --max-af 0.98  > ./pyseer_output/results/bact/bact_phenotype_ALL.tsv


pyseer --phenotypes ./pyseer_output/original_tables/bact/bact_phenotype_no_biofilm.txt --lmm \
 --pres gene_presence_absence.Rtab --similarity ./pyseer_output/phylogeny_K.tsv \
 --save-m mash_mds --max-dimensions 3 --cpu 8 --min-af 0.02 --max-af 0.98  > ./pyseer_output/results/bact/bact_phenotype_no_biofilm.tsv


pyseer --phenotypes ./pyseer_output/original_tables/bact/bact_phenotype_normal2superbio.txt --lmm \
 --pres gene_presence_absence.Rtab --similarity ./pyseer_output/phylogeny_K.tsv \
 --save-m mash_mds --max-dimensions 3 --cpu 8 --min-af 0.02 --max-af 0.98  > ./pyseer_output/results/bact/bact_phenotype_normal2superbio.tsv



### Bacterial biofilm phenotype

pyseer --phenotypes ./pyseer_output/original_tables/PG_bacterial_normal2superbio.txt --lmm \
 --pres gene_presence_absence.Rtab --similarity ./pyseer_output/phylogeny_K.tsv \
 --save-m mash_mds --max-dimensions 3 --cpu 8 --min-af 0.02 --max-af 0.98  > ./pyseer_output/results/PG_bacterial_normal2superbio.tsv






### working with unitigs


pyseer --phenotypes ./pyseer_output/original_tables/worm_phenotype_ALL_no_biofilm.txt --lmm \
 --kmers ./pyseer/unitigs/pg.unitigs.pyseer.gz --similarity ./pyseer_output/phylogeny_K.tsv \
 --save-m mash_mds --max-dimensions 3 --cpu 8 --min-af 0.02 --max-af 0.98  > ./pyseer_output/results/worm_phenotype_ALL_no_biofilm_UNITIGS.tsv





cat <(echo "#CHR SNP BP minLOG10(P) log10(p) r^2") <(paste <(sed '1d' penicillin_SNPs.txt | cut -d "_" -f 2) <(sed '1d' penicillin_SNPs.txt | cut -f 4) | awk '{p = -log($2)/log(10); print "26",".",$1,p,p,"0"}' ) | tr ' ' '\t' > penicillin_snps.plot



