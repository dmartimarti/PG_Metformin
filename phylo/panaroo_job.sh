#PBS -l walltime=48:0:0
#PBS -l select=1:ncpus=32:mem=50gb

# Load modules for any applications

module load anaconda3/personal
source activate assembly

WORK=$HOME/pangenome_study/annotations/good_quality
mkdir $WORK/annot_prokka_v1

cp $WORK/*.gff $TMPDIR

panaroo -i *.gff -o panaroo_results_1st_2nd_batch --clean-mode strict \
	-t 32 -a core --remove-invalid-genes

cp -r panaroo_results_1st_2nd_batch $WORK/
