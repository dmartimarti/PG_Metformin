#PBS -l walltime=48:0:0
#PBS -l select=1:ncpus=32:mem=24gb

# Load modules for any applications

module load anaconda3/personal
source activate panaroo

WORK=$HOME/pangenome_study/annotations/good_quality

cp $WORK/*.gff $TMPDIR

panaroo -i *.gff -o panaroo_results_1st_2nd_batch --clean-mode strict \
	-t 32 -a core --remove-invalid-genes

cp -r panaroo_results_1st_2nd_batch $WORK/
