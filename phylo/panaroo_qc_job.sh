#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=16:mem=128gb

# Load modules for any applications

module load anaconda3/personal
source activate panaroo

WORK=$HOME/pangenome_study/annotations/good_quality

mkdir results_qc

cp $WORK/*.gff $TMPDIR
cp $WORK/refseq.genomes.k21s1000.msh $TMPDIR

panaroo-qc -i *.gff -o results_qc -t 16 --graph_type all --ref_db refseq.genomes.k21s1000.msh

cp -r results_qc $WORK/