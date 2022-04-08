# parallel version of prodigal calls
# this loop takes the benefits of GNU parallel and uses all my CPUs
# within a few minutes I have the protein files from almost 800 genome assemblies

for seq in $(ls *.fasta)
do
  echo "Processing assembly $seq"
  sem -j +0 prodigal -i $seq -o $seq.genes -a $seq.faa -q
done
sem --wait