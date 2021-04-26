#PBS -l walltime=24:0:0
#PBS -l select=1:ncpus=32:mem=100gb

# Load modules for any applications

module load anaconda3/personal
source activate assembly

WORK=$HOME/Pangenomic_genomes/sequencing/raw/210401_NB501045_0331_AHYK7HBGXH/Unaligned/trimm_seqs/failed
mkdir $WORK/spades_assemblies

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

	spades.py -1 $R1paired -2 $R2paired \
	-s $R1unpaired -s $R2unpaired \
	--isolate -t 32 -m 100 \
	-k 21,33,55,77,89,99,107,115,121,127 \
	--careful

	cp -r $OUTPUT $WORK/spades_assemblies/

done