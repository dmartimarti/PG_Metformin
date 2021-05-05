#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=32:mem=124gb

# Load modules for any applications

module load anaconda3/personal
source activate iqtree

WORK=$HOME/pangenome_study/annotations/good_quality/panaroo_results_1st_2nd_batch

mkdir $WORK/phy_tree

cp $WORK/core_gene_alignment.aln $TMPDIR

iqtree -s core_gene_alignment.aln -m MFP -T AUTO --threads-max 32 -v -pre core_tree

cp core_tree* $WORK/phy_tree
