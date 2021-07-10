# download sequences from HPC

scp -r dmarti14@login.hpc.imperial.ac.uk:/rds/general/user/dmarti14/home/Pangenomic_genomes/sequencing/raw/210401_NB501045_0331_AHYK7HBGXH/Unaligned/*fastq.gz ./


#### prokka 
cd /mnt/d/MRC_Postdoc/Pangenomic/mutants_analysis/references

conda activate prokka

prokka --outdir NT12045 --genus Escherichia --prefix NT12045 --cpus 8 NT12045.fasta
prokka --outdir NT12056 --genus Escherichia --prefix NT12056 --cpus 8 NT12056.fasta
prokka --outdir NT12062 --genus Escherichia --prefix NT12062 --cpus 8 NT12062.fasta
prokka --outdir NT12321 --genus Escherichia --prefix NT12321 --cpus 8 NT12321.fasta




### trimmomatic
#example
conda activate phylo
cd /mnt/d/MRC_Postdoc/Pangenomic/mutants_analysis/raw_sequences/NT12010

for R1 in *R1*
do   	
   R2=${R1//R1.fastq.gz/R2.fastq.gz}

   R1paired=${R1//.fastq.gz/_paired.fastq.gz}
   R1unpaired=${R1//.fastq.gz/_unpaired.fastq.gz}
   R2paired=${R2//.fastq.gz/_paired.fastq.gz}
   R2unpaired=${R2//.fastq.gz/_unpaired.fastq.gz}
   
   trimmomatic PE -threads 8 -phred33 $R1 $R2 $R1paired $R1unpaired $R2paired $R2unpaired \
   ILLUMINACLIP:TruSeq3-PE.fa:4:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:140

done



## Snippy

conda activate phylo

cd /mnt/d/MRC_Postdoc/Pangenomic/mutants_analysis

# run analyses

# NT12010
snippy --outdir ./results/NT12010/NT12010_C1 --ref ./references/NT12010/NT12010.gbk --R1 ./raw_sequences/NT12010/NT12010_C1_R1_paired.fastq.gz --R2 ./raw_sequences/NT12010/NT12010_C1_R2_paired.fastq.gz --cpus 8 --report
snippy --outdir ./results/NT12010/NT12010_C2 --ref ./references/NT12010/NT12010.gbk --R1 ./raw_sequences/NT12010/NT12010_C2_R1_paired.fastq.gz --R2 ./raw_sequences/NT12010/NT12010_C2_R2_paired.fastq.gz --cpus 8 --report
snippy --outdir ./results/NT12010/NT12010_C3 --ref ./references/NT12010/NT12010.gbk --R1 ./raw_sequences/NT12010/NT12010_C3_R1_paired.fastq.gz --R2 ./raw_sequences/NT12010/NT12010_C3_R2_paired.fastq.gz --cpus 8 --report
snippy --outdir ./results/NT12010/NT12010_C4 --ref ./references/NT12010/NT12010.gbk --R1 ./raw_sequences/NT12010/NT12010_C4_R1_paired.fastq.gz --R2 ./raw_sequences/NT12010/NT12010_C4_R2_paired.fastq.gz --cpus 8 --report
snippy --outdir ./results/NT12010/NT12010_C5 --ref ./references/NT12010/NT12010.gbk --R1 ./raw_sequences/NT12010/NT12010_C5_R1_paired.fastq.gz --R2 ./raw_sequences/NT12010/NT12010_C5_R2_paired.fastq.gz --cpus 8 --report

# with input tab
snippy-multi ./raw_sequences/NT12010/input.tab  --ref ./references/NT12010/NT12010.gbk --cpus 8 > NT12010_multi.sh

snippy-core --ref 'NT12010_C1/ref.fa' NT12010_C1 NT12010_C2 NT12010_C3 NT12010_C4 NT12010_C5
mkdir NT12010_core
mv core.* ./NT12010_core

# NT12015
snippy --outdir ./results/NT12015_C1 --ref ./references/NT12015/NT12015.gbk --R1 ./raw_sequences/NT12015/NT12015_C1_R1.fastq.gz --R2 ./raw_sequences/NT12015/NT12015_C1_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/NT12015_C2 --ref ./references/NT12015/NT12015.gbk --R1 ./raw_sequences/NT12015/NT12015_C2_R1.fastq.gz --R2 ./raw_sequences/NT12015/NT12015_C2_R2.fastq.gz --cpus 8 --report
# snippy --outdir ./results/NT12015_C3 --ref ./references/NT12015/NT12015.gbk --R1 ./raw_sequences/NT12015/NT12015_C3_R1.fastq.gz --R2 ./raw_sequences/NT12015/NT12015_C3_R2.fastq.gz --cpus 8 --report # different strain
snippy --outdir ./results/NT12015_C4 --ref ./references/NT12015/NT12015.gbk --R1 ./raw_sequences/NT12015/NT12015_C4_R1.fastq.gz --R2 ./raw_sequences/NT12015/NT12015_C4_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/NT12015_C5 --ref ./references/NT12015/NT12015.gbk --R1 ./raw_sequences/NT12015/NT12015_C5_R1.fastq.gz --R2 ./raw_sequences/NT12015/NT12015_C5_R2.fastq.gz --cpus 8 --report

# NT12045
# snippy --outdir ./results/NT12045_C1 --ref ./references/NT12045/NT12045.gbk --R1 ./raw_sequences/NT12045/NT12045_C1_R1.fastq.gz --R2 ./raw_sequences/NT12045/NT12045_C1_R2.fastq.gz --cpus 8 --report # different strain
# snippy --outdir ./results/NT12045_C2 --ref ./references/NT12045/NT12045.gbk --R1 ./raw_sequences/NT12045/NT12045_C2_R1.fastq.gz --R2 ./raw_sequences/NT12045/NT12045_C2_R2.fastq.gz --cpus 8 --report # different strain
snippy --outdir ./results/NT12045_C3 --ref ./references/NT12045/NT12045.gbk --R1 ./raw_sequences/NT12045/NT12045_C3_R1.fastq.gz --R2 ./raw_sequences/NT12045/NT12045_C3_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/NT12045_C4 --ref ./references/NT12045/NT12045.gbk --R1 ./raw_sequences/NT12045/NT12045_C4_R1.fastq.gz --R2 ./raw_sequences/NT12045/NT12045_C4_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/NT12045_C5 --ref ./references/NT12045/NT12045.gbk --R1 ./raw_sequences/NT12045/NT12045_C5_R1.fastq.gz --R2 ./raw_sequences/NT12045/NT12045_C5_R2.fastq.gz --cpus 8 --report


# NT12056
# snippy --outdir ./results/NT12056_C3 --ref ./references/NT12056/NT12056.gbk --R1 ./raw_sequences/NT12056/NT12056_C3_R1.fastq.gz --R2 ./raw_sequences/NT12056/NT12056_C3_R2.fastq.gz --cpus 8 --report # different strain
snippy --outdir ./results/NT12056_C6 --ref ./references/NT12056/NT12056.gbk --R1 ./raw_sequences/NT12056/NT12056_C6_R1.fastq.gz --R2 ./raw_sequences/NT12056/NT12056_C6_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/NT12056_C5 --ref ./references/NT12056/NT12056.gbk --R1 ./raw_sequences/NT12056/NT12056_C5_R1.fastq.gz --R2 ./raw_sequences/NT12056/NT12056_C5_R2.fastq.gz --cpus 8 --report


# NT12062
snippy --outdir ./results/NT12062_C2 --ref ./references/NT12062/NT12062.gbk --R1 ./raw_sequences/NT12062/NT12062_C2_R1.fastq.gz --R2 ./raw_sequences/NT12062/NT12062_C2_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/NT12062_C3 --ref ./references/NT12062/NT12062.gbk --R1 ./raw_sequences/NT12062/NT12062_C3_R1.fastq.gz --R2 ./raw_sequences/NT12062/NT12062_C3_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/NT12062_C7 --ref ./references/NT12062/NT12062.gbk --R1 ./raw_sequences/NT12062/NT12062_C7_R1.fastq.gz --R2 ./raw_sequences/NT12062/NT12062_C7_R2.fastq.gz --cpus 8 --report

# NT12321
snippy --outdir ./results/NT12321_C4 --ref ./references/NT12321/NT12321.gbk --R1 ./raw_sequences/NT12321/NT12321_C4_R1.fastq.gz --R2 ./raw_sequences/NT12321/NT12321_C4_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/NT12321_C5 --ref ./references/NT12321/NT12321.gbk --R1 ./raw_sequences/NT12321/NT12321_C5_R1.fastq.gz --R2 ./raw_sequences/NT12321/NT12321_C5_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/NT12321_C6 --ref ./references/NT12321/NT12321.gbk --R1 ./raw_sequences/NT12321/NT12321_C6_R1.fastq.gz --R2 ./raw_sequences/NT12321/NT12321_C6_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/NT12321_C8 --ref ./references/NT12321/NT12321.gbk --R1 ./raw_sequences/NT12321/NT12321_C8_R1.fastq.gz --R2 ./raw_sequences/NT12321/NT12321_C8_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/NT12321_C27 --ref ./references/NT12321/NT12321.gbk --R1 ./raw_sequences/NT12321/NT12321_C27_R1.fastq.gz --R2 ./raw_sequences/NT12321/NT12321_C27_R2.fastq.gz --cpus 8 --report


# nissle
snippy --outdir ./results/Nissle/Nissle_P11_2_P13 --ref ./references/Nissle_Tn_parental/Nissle_Tn_parental.gbk --R1 ./raw_sequences/Nissle/Nissle_P11_2_P13_R1.fastq.gz --R2 ./raw_sequences/Nissle/Nissle_P11_2_P13_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/Nissle/Nissle_P14_2_O16 --ref ./references/Nissle_Tn_parental/Nissle_Tn_parental.gbk --R1 ./raw_sequences/Nissle/Nissle_P14_2_O16_R1.fastq.gz --R2 ./raw_sequences/Nissle/Nissle_P14_2_O16_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/Nissle/Nissle_P2_1_M17  --ref ./references/Nissle_Tn_parental/Nissle_Tn_parental.gbk --R1 ./raw_sequences/Nissle/Nissle_P2_1_M17_R1.fastq.gz --R2 ./raw_sequences/Nissle/Nissle_P2_1_M17_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/Nissle/Nissle_P3_1_G20  --ref ./references/Nissle_Tn_parental/Nissle_Tn_parental.gbk --R1 ./raw_sequences/Nissle/Nissle_P3_1_G20_R1.fastq.gz --R2 ./raw_sequences/Nissle/Nissle_P3_1_G20_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/Nissle/Nissle_P3_1_K2   --ref ./references/Nissle_Tn_parental/Nissle_Tn_parental.gbk --R1 ./raw_sequences/Nissle/Nissle_P3_1_K2_R1.fastq.gz --R2 ./raw_sequences/Nissle/Nissle_P3_1_K2_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/Nissle/Nissle_P4_1_G16  --ref ./references/Nissle_Tn_parental/Nissle_Tn_parental.gbk --R1 ./raw_sequences/Nissle/Nissle_P4_1_G16_R1.fastq.gz --R2 ./raw_sequences/Nissle/Nissle_P4_1_G16_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/Nissle/Nissle_P6_2_C1   --ref ./references/Nissle_Tn_parental/Nissle_Tn_parental.gbk --R1 ./raw_sequences/Nissle/Nissle_P6_2_C1_R1.fastq.gz --R2 ./raw_sequences/Nissle/Nissle_P6_2_C1_R2.fastq.gz --cpus 8 --report #
snippy --outdir ./results/Nissle/Nissle_P8_1_G22  --ref ./references/Nissle_Tn_parental/Nissle_Tn_parental.gbk --R1 ./raw_sequences/Nissle/Nissle_P8_1_G22_R1.fastq.gz --R2 ./raw_sequences/Nissle/Nissle_P8_1_G22_R2.fastq.gz --cpus 8 --report


