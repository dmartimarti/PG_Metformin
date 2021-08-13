# download sequences from HPC

scp -r dmarti14@login.hpc.imperial.ac.uk:/rds/general/user/dmarti14/home/Pangenomic_genomes/sequencing/raw/210401_NB501045_0331_AHYK7HBGXH/Unaligned/*fastq.gz ./


#### prokka 
cd /mnt/d/MRC_Postdoc/Pangenomic/mutants_analysis/references

conda activate prokka

prokka --outdir NT12010 --genus Escherichia --species coli --proteins /home/dani/anaconda3/envs/prokka/db/genebank/E_coli_MG1655_complete.gb --prefix NT12010 --cpus 8 NT12010.fasta
prokka --outdir NT12015 --genus Escherichia --species coli --proteins /home/dani/anaconda3/envs/prokka/db/genebank/E_coli_MG1655_complete.gb --prefix NT12015 --cpus 8 NT12015.fasta
prokka --outdir NT12045 --genus Escherichia --species coli --proteins /home/dani/anaconda3/envs/prokka/db/genebank/E_coli_MG1655_complete.gb --prefix NT12045 --cpus 8 NT12045.fasta
prokka --outdir NT12056 --genus Escherichia --species coli --proteins /home/dani/anaconda3/envs/prokka/db/genebank/E_coli_MG1655_complete.gb --prefix NT12056 --cpus 8 NT12056.fasta
prokka --outdir NT12062 --genus Escherichia --species coli --proteins /home/dani/anaconda3/envs/prokka/db/genebank/E_coli_MG1655_complete.gb --prefix NT12062 --cpus 8 NT12062.fasta
prokka --outdir NT12321 --genus Escherichia --species coli --proteins /home/dani/anaconda3/envs/prokka/db/genebank/E_coli_MG1655_complete.gb --prefix NT12321 --cpus 8 NT12321.fasta
prokka --outdir OP50 --genus Escherichia --species coli --proteins /home/dani/anaconda3/envs/prokka/db/genebank/E_coli_MG1655_complete.gb --prefix OP50 --cpus 8 OP50.fasta
prokka --outdir Nissle_Tn_parental --genus Escherichia --species coli --usegenus --proteins /home/dani/anaconda3/envs/prokka/db/genebank/E_coli_MG1655_complete.gb --prefix Nissle_Tn_parental --cpus 8 Nissle_Tn_parental.fasta

# for Nissle tn seq mutants

cd /mnt/d/MRC_Postdoc/Pangenomic/mutants_analysis/references/Nissle_references

prokka --outdir Nissle_Tn_parental --genus Escherichia --species coli --proteins /home/dani/anaconda3/envs/prokka/db/genebank/E_coli_Nissle_1917_ASM71459v1.gbff \
 --prefix Nissle_Tn_parental --cpus 8 Nissle_Tn_parental.fasta

prokka --outdir Nissle_P11_2_P13 --genus Escherichia --species coli --proteins /home/dani/anaconda3/envs/prokka/db/genebank/E_coli_Nissle_1917_ASM71459v1.gbff \
 --prefix Nissle_P11_2_P13 --cpus 8 Nissle_P11_2_P13.fasta

prokka --outdir Nissle_P2_1_M17 --genus Escherichia --species coli --proteins /home/dani/anaconda3/envs/prokka/db/genebank/E_coli_Nissle_1917_ASM71459v1.gbff \
 --prefix Nissle_P2_1_M17 --cpus 8 Nissle_P2_1_M17.fasta

prokka --outdir Nissle_P3_1_G20 --genus Escherichia --species coli --proteins /home/dani/anaconda3/envs/prokka/db/genebank/E_coli_Nissle_1917_ASM71459v1.gbff \
 --prefix Nissle_P3_1_G20 --cpus 8 Nissle_P3_1_G20.fasta

prokka --outdir Nissle_P3_1_K2 --genus Escherichia --species coli --proteins /home/dani/anaconda3/envs/prokka/db/genebank/E_coli_Nissle_1917_ASM71459v1.gbff \
 --prefix Nissle_P3_1_K2 --cpus 8 Nissle_P3_1_K2.fasta

prokka --outdir Nissle_P4_1_G16 --genus Escherichia --species coli --proteins /home/dani/anaconda3/envs/prokka/db/genebank/E_coli_Nissle_1917_ASM71459v1.gbff \
 --prefix Nissle_P4_1_G16 --cpus 8 Nissle_P4_1_G16.fasta

prokka --outdir Nissle_P6_2_C1 --genus Escherichia --species coli --proteins /home/dani/anaconda3/envs/prokka/db/genebank/E_coli_Nissle_1917_ASM71459v1.gbff \
 --prefix Nissle_P6_2_C1 --cpus 8 Nissle_P6_2_C1.fasta

prokka --outdir Nissle_P8_1_G22 --genus Escherichia --species coli --proteins /home/dani/anaconda3/envs/prokka/db/genebank/E_coli_Nissle_1917_ASM71459v1.gbff \
 --prefix Nissle_P8_1_G22 --cpus 8 Nissle_P8_1_G22.fasta

prokka --outdir Nissle_P14_2_O16 --genus Escherichia --species coli --proteins /home/dani/anaconda3/envs/prokka/db/genebank/E_coli_Nissle_1917_ASM71459v1.gbff \
 --prefix Nissle_P14_2_O16 --cpus 8 Nissle_P14_2_O16.fasta



### trimmomatic
#example
conda activate phylo
cd /mnt/d/MRC_Postdoc/Pangenomic/mutants_analysis/raw_sequences/Nissle

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

mkdir filtered

mv *paired* ./filtered



## Snippy

conda activate phylo

cd /mnt/d/MRC_Postdoc/Pangenomic/mutants_analysis

# run analyses

# NT12010
snippy --outdir ./results/NT12010/NT12010_C1 --ref ./references/NT12010/NT12010.gbk --R1 ./raw_sequences/NT12010/NT12010_C1_R1.fastq.gz --R2 ./raw_sequences/NT12010/NT12010_C1_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/NT12010/NT12010_C2 --ref ./references/NT12010/NT12010.gbk --R1 ./raw_sequences/NT12010/NT12010_C2_R1.fastq.gz --R2 ./raw_sequences/NT12010/NT12010_C2_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/NT12010/NT12010_C3 --ref ./references/NT12010/NT12010.gbk --R1 ./raw_sequences/NT12010/NT12010_C3_R1.fastq.gz --R2 ./raw_sequences/NT12010/NT12010_C3_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/NT12010/NT12010_C4 --ref ./references/NT12010/NT12010.gbk --R1 ./raw_sequences/NT12010/NT12010_C4_R1.fastq.gz --R2 ./raw_sequences/NT12010/NT12010_C4_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/NT12010/NT12010_C5 --ref ./references/NT12010/NT12010.gbk --R1 ./raw_sequences/NT12010/NT12010_C5_R1.fastq.gz --R2 ./raw_sequences/NT12010/NT12010_C5_R2.fastq.gz --cpus 8 --report

# core 
cd ./results/NT12010
snippy-core --ref 'NT12010_C1/ref.fa' NT12010_C1 NT12010_C2 NT12010_C3 NT12010_C4 NT12010_C5
mkdir NT12010_core
mv core.* ./NT12010_core
cd ../../



# NT12015
snippy --outdir ./results/NT12015/NT12015_C1 --ref ./references/NT12015/NT12015.gbk --R1 ./raw_sequences/NT12015/NT12015_C1_R1.fastq.gz --R2 ./raw_sequences/NT12015/NT12015_C1_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/NT12015/NT12015_C2 --ref ./references/NT12015/NT12015.gbk --R1 ./raw_sequences/NT12015/NT12015_C2_R1.fastq.gz --R2 ./raw_sequences/NT12015/NT12015_C2_R2.fastq.gz --cpus 8 --report
# snippy --outdir ./results/NT12015/NT12015_C3 --ref ./references/NT12015/NT12015.gbk --R1 ./raw_sequences/NT12015/NT12015_C3_R1.fastq.gz --R2 ./raw_sequences/NT12015/NT12015_C3_R2.fastq.gz --cpus 8 --report # different strain
snippy --outdir ./results/NT12015/NT12015_C4 --ref ./references/NT12015/NT12015.gbk --R1 ./raw_sequences/NT12015/NT12015_C4_R1.fastq.gz --R2 ./raw_sequences/NT12015/NT12015_C4_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/NT12015/NT12015_C5 --ref ./references/NT12015/NT12015.gbk --R1 ./raw_sequences/NT12015/NT12015_C5_R1.fastq.gz --R2 ./raw_sequences/NT12015/NT12015_C5_R2.fastq.gz --cpus 8 --report

# NT12045
# snippy --outdir ./results/NT12045_C1 --ref ./references/NT12045/NT12045.gbk --R1 ./raw_sequences/NT12045/NT12045_C1_R1.fastq.gz --R2 ./raw_sequences/NT12045/NT12045_C1_R2.fastq.gz --cpus 8 --report # different strain
# snippy --outdir ./results/NT12045_C2 --ref ./references/NT12045/NT12045.gbk --R1 ./raw_sequences/NT12045/NT12045_C2_R1.fastq.gz --R2 ./raw_sequences/NT12045/NT12045_C2_R2.fastq.gz --cpus 8 --report # different strain
snippy --outdir ./results/NT12045/NT12045_C3 --ref ./references/NT12045/NT12045.gbk --R1 ./raw_sequences/NT12045/NT12045_C3_R1.fastq.gz --R2 ./raw_sequences/NT12045/NT12045_C3_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/NT12045/NT12045_C4 --ref ./references/NT12045/NT12045.gbk --R1 ./raw_sequences/NT12045/NT12045_C4_R1.fastq.gz --R2 ./raw_sequences/NT12045/NT12045_C4_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/NT12045/NT12045_C5 --ref ./references/NT12045/NT12045.gbk --R1 ./raw_sequences/NT12045/NT12045_C5_R1.fastq.gz --R2 ./raw_sequences/NT12045/NT12045_C5_R2.fastq.gz --cpus 8 --report


# NT12056
# snippy --outdir ./results/NT12056_C3 --ref ./references/NT12056/NT12056.gbk --R1 ./raw_sequences/NT12056/NT12056_C3_R1.fastq.gz --R2 ./raw_sequences/NT12056/NT12056_C3_R2.fastq.gz --cpus 8 --report # different strain
snippy --outdir ./results/NT12056/NT12056_C6 --ref ./references/NT12056/NT12056.gbk --R1 ./raw_sequences/NT12056/NT12056_C6_R1.fastq.gz --R2 ./raw_sequences/NT12056/NT12056_C6_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/NT12056/NT12056_C5 --ref ./references/NT12056/NT12056.gbk --R1 ./raw_sequences/NT12056/NT12056_C5_R1.fastq.gz --R2 ./raw_sequences/NT12056/NT12056_C5_R2.fastq.gz --cpus 8 --report


# NT12062
snippy --outdir ./results/NT12062/NT12062_C2 --ref ./references/NT12062/NT12062.gbk --R1 ./raw_sequences/NT12062/NT12062_C2_R1.fastq.gz --R2 ./raw_sequences/NT12062/NT12062_C2_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/NT12062/NT12062_C3 --ref ./references/NT12062/NT12062.gbk --R1 ./raw_sequences/NT12062/NT12062_C3_R1.fastq.gz --R2 ./raw_sequences/NT12062/NT12062_C3_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/NT12062/NT12062_C7 --ref ./references/NT12062/NT12062.gbk --R1 ./raw_sequences/NT12062/NT12062_C7_R1.fastq.gz --R2 ./raw_sequences/NT12062/NT12062_C7_R2.fastq.gz --cpus 8 --report

# NT12321
snippy --outdir ./results/NT12321/NT12321_C4 --ref ./references/NT12321/NT12321.gbk --R1 ./raw_sequences/NT12321/NT12321_C4_R1.fastq.gz --R2 ./raw_sequences/NT12321/NT12321_C4_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/NT12321/NT12321_C5 --ref ./references/NT12321/NT12321.gbk --R1 ./raw_sequences/NT12321/NT12321_C5_R1.fastq.gz --R2 ./raw_sequences/NT12321/NT12321_C5_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/NT12321/NT12321_C6 --ref ./references/NT12321/NT12321.gbk --R1 ./raw_sequences/NT12321/NT12321_C6_R1.fastq.gz --R2 ./raw_sequences/NT12321/NT12321_C6_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/NT12321/NT12321_C8 --ref ./references/NT12321/NT12321.gbk --R1 ./raw_sequences/NT12321/NT12321_C8_R1.fastq.gz --R2 ./raw_sequences/NT12321/NT12321_C8_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/NT12321/NT12321_C27 --ref ./references/NT12321/NT12321.gbk --R1 ./raw_sequences/NT12321/NT12321_C27_R1.fastq.gz --R2 ./raw_sequences/NT12321/NT12321_C27_R2.fastq.gz --cpus 8 --report


# nissle
snippy --outdir ./results/Nissle/Nissle_P11_2_P13 --ref ./references/Nissle_Tn_parental/Nissle_Tn_parental.gbk --R1 ./raw_sequences/Nissle/Nissle_P11_2_P13_R1.fastq.gz --R2 ./raw_sequences/Nissle/Nissle_P11_2_P13_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/Nissle/Nissle_P14_2_O16 --ref ./references/Nissle_Tn_parental/Nissle_Tn_parental.gbk --R1 ./raw_sequences/Nissle/Nissle_P14_2_O16_R1.fastq.gz --R2 ./raw_sequences/Nissle/Nissle_P14_2_O16_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/Nissle/Nissle_P2_1_M17  --ref ./references/Nissle_Tn_parental/Nissle_Tn_parental.gbk --R1 ./raw_sequences/Nissle/Nissle_P2_1_M17_R1.fastq.gz --R2 ./raw_sequences/Nissle/Nissle_P2_1_M17_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/Nissle/Nissle_P3_1_G20  --ref ./references/Nissle_Tn_parental/Nissle_Tn_parental.gbk --R1 ./raw_sequences/Nissle/Nissle_P3_1_G20_R1.fastq.gz --R2 ./raw_sequences/Nissle/Nissle_P3_1_G20_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/Nissle/Nissle_P3_1_K2   --ref ./references/Nissle_Tn_parental/Nissle_Tn_parental.gbk --R1 ./raw_sequences/Nissle/Nissle_P3_1_K2_R1.fastq.gz --R2 ./raw_sequences/Nissle/Nissle_P3_1_K2_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/Nissle/Nissle_P4_1_G16  --ref ./references/Nissle_Tn_parental/Nissle_Tn_parental.gbk --R1 ./raw_sequences/Nissle/Nissle_P4_1_G16_R1.fastq.gz --R2 ./raw_sequences/Nissle/Nissle_P4_1_G16_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/Nissle/Nissle_P6_2_C1   --ref ./references/Nissle_Tn_parental/Nissle_Tn_parental.gbk --R1 ./raw_sequences/Nissle/Nissle_P6_2_C1_R1.fastq.gz --R2 ./raw_sequences/Nissle/Nissle_P6_2_C1_R2.fastq.gz --cpus 8 --report #
snippy --outdir ./results/Nissle/Nissle_P8_1_G22  --ref ./references/Nissle_Tn_parental/Nissle_Tn_parental.gbk --R1 ./raw_sequences/Nissle/Nissle_P8_1_G22_R1.fastq.gz --R2 ./raw_sequences/Nissle/Nissle_P8_1_G22_R2.fastq.gz --cpus 8 --report


# M strains

snippy --outdir ./results/M_strains/M1 --ref ./references/OP50/OP50.gbk --R1 ./raw_sequences/M_strains/M1_R1.fastq.gz --R2 ./raw_sequences/M_strains/M1_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/M_strains/M2 --ref ./references/OP50/OP50.gbk --R1 ./raw_sequences/M_strains/M2_R1.fastq.gz --R2 ./raw_sequences/M_strains/M2_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/M_strains/M3 --ref ./references/OP50/OP50.gbk --R1 ./raw_sequences/M_strains/M3_R1.fastq.gz --R2 ./raw_sequences/M_strains/M3_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/M_strains/M4 --ref ./references/OP50/OP50.gbk --R1 ./raw_sequences/M_strains/M4_R1.fastq.gz --R2 ./raw_sequences/M_strains/M4_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/M_strains/M5 --ref ./references/OP50/OP50.gbk --R1 ./raw_sequences/M_strains/M5_R1.fastq.gz --R2 ./raw_sequences/M_strains/M5_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/M_strains/M6 --ref ./references/OP50/OP50.gbk --R1 ./raw_sequences/M_strains/M6_R1.fastq.gz --R2 ./raw_sequences/M_strains/M6_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/M_strains/M7 --ref ./references/OP50/OP50.gbk --R1 ./raw_sequences/M_strains/M7_R1.fastq.gz --R2 ./raw_sequences/M_strains/M7_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/M_strains/M8 --ref ./references/OP50/OP50.gbk --R1 ./raw_sequences/M_strains/M8_R1.fastq.gz --R2 ./raw_sequences/M_strains/M8_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/M_strains/M9 --ref ./references/OP50/OP50.gbk --R1 ./raw_sequences/M_strains/M9_R1.fastq.gz --R2 ./raw_sequences/M_strains/M9_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/M_strains/M10 --ref ./references/OP50/OP50.gbk --R1 ./raw_sequences/M_strains/M10_R1.fastq.gz --R2 ./raw_sequences/M_strains/M10_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/M_strains/M11 --ref ./references/OP50/OP50.gbk --R1 ./raw_sequences/M_strains/M11_R1.fastq.gz --R2 ./raw_sequences/M_strains/M11_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/M_strains/M12 --ref ./references/OP50/OP50.gbk --R1 ./raw_sequences/M_strains/M12_R1.fastq.gz --R2 ./raw_sequences/M_strains/M12_R2.fastq.gz --cpus 8 --report

# P strains
snippy --outdir ./results/P_strains/P1 --ref ./references/OP50/OP50.gbk --R1 ./raw_sequences/P_strains/P1_R1.fastq.gz --R2 ./raw_sequences/P_strains/P1_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/P_strains/P2 --ref ./references/OP50/OP50.gbk --R1 ./raw_sequences/P_strains/P2_R1.fastq.gz --R2 ./raw_sequences/P_strains/P2_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/P_strains/P3 --ref ./references/OP50/OP50.gbk --R1 ./raw_sequences/P_strains/P3_R1.fastq.gz --R2 ./raw_sequences/P_strains/P3_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/P_strains/P4 --ref ./references/OP50/OP50.gbk --R1 ./raw_sequences/P_strains/P4_R1.fastq.gz --R2 ./raw_sequences/P_strains/P4_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/P_strains/P5 --ref ./references/OP50/OP50.gbk --R1 ./raw_sequences/P_strains/P5_R1.fastq.gz --R2 ./raw_sequences/P_strains/P5_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/P_strains/P6 --ref ./references/OP50/OP50.gbk --R1 ./raw_sequences/P_strains/P6_R1.fastq.gz --R2 ./raw_sequences/P_strains/P6_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/P_strains/P7 --ref ./references/OP50/OP50.gbk --R1 ./raw_sequences/P_strains/P7_R1.fastq.gz --R2 ./raw_sequences/P_strains/P7_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/P_strains/P8 --ref ./references/OP50/OP50.gbk --R1 ./raw_sequences/P_strains/P8_R1.fastq.gz --R2 ./raw_sequences/P_strains/P8_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/P_strains/P9 --ref ./references/OP50/OP50.gbk --R1 ./raw_sequences/P_strains/P9_R1.fastq.gz --R2 ./raw_sequences/P_strains/P9_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/P_strains/P10 --ref ./references/OP50/OP50.gbk --R1 ./raw_sequences/P_strains/P10_R1.fastq.gz --R2 ./raw_sequences/P_strains/P10_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/P_strains/P11 --ref ./references/OP50/OP50.gbk --R1 ./raw_sequences/P_strains/P11_R1.fastq.gz --R2 ./raw_sequences/P_strains/P11_R2.fastq.gz --cpus 11 --report
snippy --outdir ./results/P_strains/P12 --ref ./references/OP50/OP50.gbk --R1 ./raw_sequences/P_strains/P12_R1.fastq.gz --R2 ./raw_sequences/P_strains/P12_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/P_strains/P13 --ref ./references/OP50/OP50.gbk --R1 ./raw_sequences/P_strains/P13_R1.fastq.gz --R2 ./raw_sequences/P_strains/P13_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/P_strains/P14 --ref ./references/OP50/OP50.gbk --R1 ./raw_sequences/P_strains/P14_R1.fastq.gz --R2 ./raw_sequences/P_strains/P14_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/P_strains/P15 --ref ./references/OP50/OP50.gbk --R1 ./raw_sequences/P_strains/P15_R1.fastq.gz --R2 ./raw_sequences/P_strains/P15_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/P_strains/P16 --ref ./references/OP50/OP50.gbk --R1 ./raw_sequences/P_strains/P16_R1.fastq.gz --R2 ./raw_sequences/P_strains/P16_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/P_strains/P17 --ref ./references/OP50/OP50.gbk --R1 ./raw_sequences/P_strains/P17_R1.fastq.gz --R2 ./raw_sequences/P_strains/P17_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/P_strains/P18 --ref ./references/OP50/OP50.gbk --R1 ./raw_sequences/P_strains/P18_R1.fastq.gz --R2 ./raw_sequences/P_strains/P18_R2.fastq.gz --cpus 8 --report
snippy --outdir ./results/P_strains/P19 --ref ./references/OP50/OP50.gbk --R1 ./raw_sequences/P_strains/P19_R1.fastq.gz --R2 ./raw_sequences/P_strains/P19_R2.fastq.gz --cpus 8 --report




### JBrowse 2

cd /mnt/d/MRC_Postdoc/Pangenomic/mutants_analysis/results/P_strains/P16

jbrowse create jbrowse2
cd jbrowse2

jbrowse add-assembly ../reference/ref.fa --load copy 
jbrowse add-track ../snps.bam --load copy
tabix ../snps.vcf.gz
jbrowse add-track ../snps.vcf.gz --load copy

gt gff3 -sortlines -tidy -retainids ../reference/ref.gff > ../reference/ref.sorted.gff
bgzip ../reference/ref.sorted.gff
tabix ../reference/ref.sorted.gff.gz
jbrowse add-track ../reference/ref.sorted.gff.gz --load copy
# run jbrowse 2
npx serve .









### Genome assembly of Nissle with Nissle reference

cd /mnt/d/MRC_Postdoc/Pangenomic/mutants_analysis/raw_sequences/Nissle/filtered

# ref = GCF_000714595.1_ASM71459v1_genomic.fna

conda activate assembly

spades.py -1 Nissle_Tn_parental_R1_001_paired.fastq.gz -2 Nissle_Tn_parental_R2_001_paired.fastq.gz \
	-s Nissle_Tn_parental_R1_001_unpaired.fastq.gz -s Nissle_Tn_parental_R2_001_unpaired.fastq.gz \
	--isolate -t 8 -m 16 --trusted-contigs ASM71459v1.fasta \
	-k 21,33,55,77,89,99,107,115,121,127 -o Nissle_parental 

