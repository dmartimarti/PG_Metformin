#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=16:mem=64gb

# Load modules for any applications

module load anaconda3/personal
source activate iqtree

WORK=$HOME/pangenome_study/complete/annotations/no_evo/panaroo_results_noEVO

mkdir $WORK/phy_fast_tree

cp $WORK/core_gene_alignment.aln $TMPDIR

iqtree -s core_gene_alignment.aln -fast -m GTR -T AUTO --threads-max 16 -v -pre core_fast_tree

cp core_fast_tree* $WORK/phy_fast_tree
