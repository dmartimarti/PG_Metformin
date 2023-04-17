## let's try to use unitigs with Pyseer and annotate the results to reference genomes


# to run pyseer, se need to go to this path and run this command 

# this works for my Macbook, as I have both miniforge (M1 native) and miniconda (rossetta version)
source ~/start_miniconda.sh 

conda activate pyseer

cd /Users/danmarti/Documents/MRC_postdoc/Pangenomic/pangenome_analysis/ALL/phylo_analysis/panaroo_results_noEVO

# pyseer command line, pay attention that we save the kmer patterns as well
pyseer --phenotypes ./pyseer_output/original_tables/worm_phenotype_ALL_no_biofilm.txt --lmm \
 --kmers ./pyseer/unitigs/pg.unitigs.pyseer.gz --similarity ./pyseer_output/phylogeny_K.tsv \
 --save-m mash_mds --max-dimensions 3 --cpu 8 --min-af 0.02 --max-af 0.98  \
 --output-patterns ./pyseer_output/results/wp_ALL_no_biofilm_kmer_patterns.txt > ./pyseer_output/results/worm_phenotype_ALL_no_biofilm_UNITIGS.tsv

# the output file is worm_phenotype_ALL_no_biofilm_UNITIGS.tsv

# calculate significance threshold
python pyseer/scripts/count_patterns.py ./pyseer_output/results/wp_ALL_no_biofilm_kmer_patterns.txt 

##### output
#Patterns:    1614173
#Threshold:  3.10E-08


# draw a QQ plot
python pyseer/scripts/qq_plot.py ./pyseer_output/results/worm_phenotype_ALL_no_biofilm_UNITIGS.tsv ./pyseer_output/results/worm_phenotype_ALL_no_biofilm_UNITIGS 


# Interpreting significant k-mers

cd /Users/danmarti/Documents/MRC_postdoc/Pangenomic/pangenome_analysis/ALL/phylo_analysis/panaroo_results_noEVO/pyseer_output/results/

cat <(head -1 worm_phenotype_ALL_no_biofilm_UNITIGS.tsv) <(awk '$4<3.10E-04 {print $0}' worm_phenotype_ALL_no_biofilm_UNITIGS.tsv) > wp_ALL_no_biofilm_significant_kmers.txt


cat <(head -1 worm_phenotype_ALL_no_biofilm_UNITIGS.tsv) <(awk '$3<3.10E-05 || $4<3.10E-05 {print $0}' worm_phenotype_ALL_no_biofilm_UNITIGS.tsv) > wp_ALL_no_biofilm_significant_kmers.txt


# Mapping to a single reference
python ../../pyseer/phandango_mapper-runner.py  wp_ALL_no_biofilm_significant_kmers.txt ../ecoli_reference/GCF_000005845.2_ASM584v2_genomic.fna ./wp_ALL_no_biofilm_significant_kmers.plot

python ../../pyseer/phandango_mapper-runner.py  worm_phenotype_ALL_no_biofilm_UNITIGS.tsv ../ecoli_reference/GCF_000005845.2_ASM584v2_genomic.fna ./wp_ALL_no_biofilm_unitigs.plot


# annotate the kmers
python  ../pyseer/annotate_hits_pyseer-runner.py results/wp_ALL_no_biofilm_significant_kmers.txt ref_complete.txt results/wp_ALL_no_biofilm_annotated_kmers.txt

# summarise 
python ../pyseer/scripts/summarise_annotations.py  results/wp_ALL_no_biofilm_annotated_kmers.txt > results/gene_hits.txt




