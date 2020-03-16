
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