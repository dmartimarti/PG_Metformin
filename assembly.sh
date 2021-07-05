
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
for folder in P5*
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


## GENOMIC DISTANCES

phylonium *.fasta > distances.txt
# REMEMBER TO REMOVE FIRST ROW FROM THE FILE




## COPY ANNOTATIONS MADE FROM PROKKA INTO MY COMPUTER
cd /mnt/d/MRC_Postdoc/Pangenomic/sequencing/Big_sequencing/unicycler_assemblies/annot_prokka

# change the folder for the one you want to download
scp -r dmarti14@login.hpc.imperial.ac.uk:/rds/general/user/dmarti14/home/Pangenomic_genomes/sequencing/raw/210401_NB501045_0331_AHYK7HBGXH/all_assemblies/annot_prokka_v1/NT12060 ./




### PHYLOGROUP ANALYSIS

../ClermonTyping/clermonTyping.sh --fasta 1.2.fasta@1.3.fasta@1.fasta@10.fasta@100.fasta@101.fasta@102.fasta@103.fasta@104.fasta@105.fasta@106.fasta@107.fasta@108.fasta@109.fasta@11.fasta@110.fasta@111.fasta@112.fasta@113.fasta@114.fasta@115.fasta@116.fasta@117.fasta@118.fasta@119.fasta@120.fasta@121.fasta@122.fasta@123.fasta@124.fasta@125.fasta@126.fasta@127.fasta@128.fasta@129.fasta@130.fasta@131.fasta@132.fasta@133.fasta@134.fasta@135.fasta@136.fasta@137.fasta@138.fasta@139.fasta@14.fasta@140.fasta@141.fasta@142.fasta@143.fasta@144.fasta@145.fasta@146.fasta@147.fasta@148.fasta@149.fasta@15.fasta@150.fasta@151.fasta@152.fasta@154.fasta@155.fasta@159.fasta@16.fasta@160.fasta@165.fasta@166.fasta@167.fasta@168.fasta@169.fasta@17.fasta@172.fasta@173.fasta@174.fasta@175.fasta@176.fasta@178.fasta@18.fasta@180.fasta@183.fasta@184.fasta@185.fasta@186.fasta@19.fasta@191.fasta@193.fasta@197.fasta@199.fasta@2.3.fasta@2.fasta@20.fasta@200.fasta@201.fasta@202.fasta@203.fasta@204.fasta@205.fasta@206.fasta@207.fasta@208.fasta@209.fasta@21.fasta@210.fasta@211.fasta@212.fasta@214.fasta@215.fasta@216.fasta@217.fasta@219.fasta@22.fasta@220.fasta@222.fasta@223.fasta@224.fasta@225.fasta@228.fasta@23.fasta@231.fasta@232.fasta@236.fasta@237.fasta@239.fasta@24.fasta@240.fasta@244.fasta@25.fasta@26.fasta@27.fasta@28.fasta@29.fasta@3.fasta@30.fasta@31.fasta@32.fasta@33.fasta@34.fasta@35.fasta@36.fasta@37.fasta@38.fasta@39.fasta@40.fasta@41.fasta@42.fasta@43.fasta@44.fasta@45.fasta@46.fasta@47.fasta@48.fasta@49.fasta@50.fasta@51.fasta@52.fasta@53.fasta@54.fasta@55.fasta@56.fasta@57.fasta@58.fasta@59.fasta@6.fasta@60.fasta@61.fasta@62.fasta@63.fasta@64.fasta@65.fasta@66.fasta@67.fasta@68.fasta@69.fasta@70.fasta@71.fasta@73.fasta@74.fasta@75.fasta@76.fasta@77.fasta@78.fasta@79.fasta@8.fasta@80.fasta@81.fasta@82.fasta@83.fasta@84.fasta@85.fasta@86.fasta@87.fasta@88.fasta@89.fasta@9.fasta@90.fasta@91.fasta@92.fasta@93.fasta@94.fasta@95.fasta@96.fasta@97.fasta@98.fasta@99.fasta@FOC_10.fasta@FOC_3.fasta@FOC_5.1.fasta@FOC_7.1.fasta@NT12010.fasta@NT12010_C1.fasta@NT12010_C2.fasta@NT12010_C3.fasta@NT12010_C4.fasta@NT12010_C5.fasta@NT12015.fasta@NT12015_C1.fasta@NT12015_C2.fasta@NT12015_C3.fasta@NT12015_C4.fasta@NT12015_C5.fasta@NT12045.fasta@NT12045_C1.fasta@NT12045_C2.fasta@NT12045_C3.fasta@NT12045_C4.fasta@NT12045_C5.fasta@NT12050.fasta@NT12056.fasta@NT12056_C3.fasta@NT12056_C5.fasta@NT12056_C6.fasta@NT12060.fasta@NT12062.fasta@NT12062_C2.fasta@NT12062_C3.fasta@NT12062_C7.fasta@NT12078.fasta@NT12129.fasta@NT12142.fasta@NT12149.fasta@NT12198.fasta@NT12200.fasta@NT12201.fasta@NT12203.fasta@NT12206.fasta@NT12207.fasta@NT12208.fasta@NT12209.fasta@NT12214.fasta@NT12216.fasta@NT12222.fasta@NT12225.fasta@NT12227.fasta@NT12228.fasta@NT12229.fasta@NT12231.fasta@NT12234.fasta@NT12235.fasta@NT12236.fasta@NT12237.fasta@NT12238.fasta@NT12239.fasta@NT12240.fasta@NT12241.fasta@NT12242.fasta@NT12243.fasta@NT12244.fasta@NT12245.fasta@NT12246.fasta@NT12247.fasta@NT12248.fasta@NT12249.fasta@NT12250.fasta@NT12251.fasta@NT12252.fasta@NT12253.fasta@NT12254.fasta@NT12255.fasta@NT12256.fasta@NT12258.fasta@NT12259.fasta@NT12260.fasta@NT12261.fasta@NT12263.fasta@NT12264.fasta@NT12265.fasta@NT12266.fasta@NT12267.fasta@NT12268.fasta@NT12269.fasta@NT12270.fasta@NT12271.fasta@NT12272.fasta@NT12273.fasta@NT12274.fasta@NT12275.fasta@NT12276.fasta@NT12277.fasta@NT12278.fasta@NT12279.fasta@NT12280.fasta@NT12281.fasta@NT12282.fasta@NT12283.fasta@NT12284.fasta@NT12287.fasta@NT12289.fasta@NT12290.fasta@NT12291.fasta@NT12293.fasta@NT12295.fasta@NT12296.fasta@NT12297.fasta@NT12300.fasta@NT12301.fasta@NT12304.fasta@NT12307.fasta@NT12312.fasta@NT12313.fasta@NT12316.fasta@NT12318.fasta@NT12321.fasta@NT12321_C27.fasta@NT12321_C4.fasta@NT12321_C5.fasta@NT12321_C6.fasta@NT12321_C8.fasta@NT12322.fasta@NT12323.fasta@NT12324.fasta@NT12339.fasta@NT12344.fasta@NT12345.fasta@NT12348.fasta@NT12349.fasta@NT12350.fasta@NT12351.fasta@NT12352.fasta@NT12354.fasta@NT12356.fasta@NT12357.fasta@NT12368.fasta@NT12369.fasta@NT12371.fasta@NT12373.fasta@NT12375.fasta@NT12376.fasta@NT12377.fasta@NT12379.fasta@NT12380.fasta@NT12381.fasta@NT12385.fasta@NT12386.fasta@NT12395.fasta@NT12397.fasta@NT12399.fasta@NT12401.fasta@NT12402.fasta@NT12405.fasta@NT12409.fasta@NT12410.fasta@NT12411.fasta@NT12415.fasta@NT12416.fasta@NT12417.fasta@NT12419.fasta@NT12428.fasta@NT12430.fasta@NT12431.fasta@NT12433.fasta@NT12437.fasta@NT12439.fasta@NT12441.fasta@NT12443.fasta@NT12446.fasta@NT12613.fasta@OP50.fasta@P1.fasta@P16.fasta@P2.fasta@P3.fasta@P4.fasta@SPC_3.1.fasta@SPC_4.fasta --threshold 250

