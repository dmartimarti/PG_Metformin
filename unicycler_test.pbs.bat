#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=32:mem=30gb

# Load modules for any applications

module load anaconda3/personal
source activate assembly

WORK=$HOME/Pangenomic_genomes/sequencing/raw/210401_NB501045_0331_AHYK7HBGXH/test/trimm_seqs_minlen130
mkdir $WORK/unicycler_batch_example

for R1paired in $WORK/*R1_001_paired*
do   	
	R2paired=${R1paired//R1_001_paired.fastq.gz/R2_001_paired.fastq.gz}
	R1unpaired=${R1paired//R1_001_paired.fastq.gz/R1_001_unpaired.fastq.gz}   	
	R2unpaired=${R1paired//R1_001_paired.fastq.gz/R2_001_unpaired.fastq.gz}
	OUTLONG=${R1paired##*/}
    OUTPUT=${OUTLONG::-23}

    cp $R1paired $TMPDIR
    cp $R2paired $TMPDIR
    cp $R1unpaired $TMPDIR
    cp $R2unpaired $TMPDIR

    unicycler -1 $R1paired -2 $R2paired \
	-s $R1unpaired -s $R2unpaired \
	-o $OUTPUT -t 32 --no_correct

	cp -r $OUTPUT $WORK/unicycler_batch_example/

done
