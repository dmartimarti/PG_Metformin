
# B1
# load path where the sequences are, and cat them

cd /Users/dmarti14/Documents/MRC_Postdoc/Projects/Pangenomic/sequencing/B1/fastq_pass 

cat *.fastq > input.fastq 

filtlong --min_length 1000 --keep_percent 90 input.fastq | gzip > output.fastq

unicycler -l output.fastq -o ../Unicycler_assembly -t 6 



# C1

cd /Users/dmarti14/Documents/MRC_Postdoc/Projects/Pangenomic/sequencing/C1/fastq_pass

cat *.fastq > input.fastq 

filtlong --min_length 500 --keep_percent 95 input.fastq | gzip > output.fastq

unicycler -l output.fastq -o ../Unicycler_assembly -t 6 




####
# These command lines are meant for the assembly of the big bulk of genomes from the experiments
# I'm using the HPC at Imperial as well as my own computer to do everything


# this code collects all genomes from a batch folder 
# remember to change P (plate) and batch (1 to 8)
# P1: batch1 and batch2
# P2: batch3 and batch4
# P3: batch5 and batch6
# P4: batch7 and batch8

cd /rds/general/user/dmarti14/home/Pangenomic_genomes/sequencing/raw/210401_NB501045_0331_AHYK7HBGXH/Unaligned/trimm_seqs/P1/batch1/unicycler_batch
mkdir assemblies
for folder in P1*
do
	cp ./$folder/assembly.fasta ./assemblies/$folder.fasta
done

# once we have all genomes in a single folder, copy them into my computer to run quast (it doesn't run in the HPC for some reason)

# scp to copy fasta files into my folder
# my machine
cd /mnt/d/MRC_Postdoc/Pangenomic/sequencing/Big_sequencing/unicycler_assemblies/P1
scp -r dmarti14@login.hpc.imperial.ac.uk:/rds/general/user/dmarti14/home/Pangenomic_genomes/sequencing/raw/210401_NB501045_0331_AHYK7HBGXH/Unaligned/trimm_seqs/P4/batch8/unicycler_batch/assemblies/*fasta ./


cd /mnt/d/MRC_Postdoc/Pangenomic/sequencing/Big_sequencing/unicycler_assemblies/unicycler_files
scp -r dmarti14@login.hpc.imperial.ac.uk:/rds/general/user/dmarti14/home/Pangenomic_genomes/sequencing/raw/210401_NB501045_0331_AHYK7HBGXH/Unaligned/trimm_seqs/P*/batch*/unicycler_batch/P* ./

# quast
conda activate assembly
quast *fasta -t 10 -o quast_quality




for folder in *
do
	cp ./$folder/*.fasta ./all_assemblies/
done