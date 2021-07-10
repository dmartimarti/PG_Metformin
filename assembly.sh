
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





### NEW ANALYSIS
# this is done to the final list of files, both for non-evo and for the folder including evo experiment strains

cd /mnt/d/MRC_Postdoc/Pangenomic/pangenome_analysis/ALL/phylo_analysis/assemblies/no_evo

../ClermonTyping/clermonTyping.sh --fasta 1.2.fasta@1.3.fasta@1.fasta@100.fasta@101.fasta@102.fasta@103.fasta@104.fasta@105.fasta@106.fasta@107.fasta@108.fasta@109.fasta@110.fasta@111.fasta@112.fasta@113.fasta@114.fasta@115.fasta@116.fasta@117.fasta@118.fasta@119.fasta@120.fasta@121.fasta@123.fasta@124.fasta@125.fasta@126.fasta@127.fasta@128.fasta@129.fasta@130.fasta@131.fasta@132.fasta@133.fasta@134.fasta@135.fasta@136.fasta@137.fasta@138.fasta@139.fasta@14.fasta@140.fasta@141.fasta@142.fasta@143.fasta@144.fasta@145.fasta@146.fasta@147.fasta@148.fasta@149.fasta@15.fasta@150.fasta@151.fasta@152.fasta@153.fasta@155.fasta@156.fasta@157.fasta@158.fasta@159.fasta@16.fasta@160.fasta@161.fasta@162.fasta@163.fasta@164.fasta@165.fasta@166.fasta@167.fasta@168.fasta@169.fasta@170.fasta@171.fasta@172.fasta@173.fasta@174.fasta@175.fasta@176.fasta@177.fasta@179.fasta@18.fasta@180.fasta@181.fasta@182.fasta@183.fasta@184.fasta@185.fasta@186.fasta@187.fasta@188.fasta@189.fasta@19.fasta@191.fasta@192.fasta@193.fasta@194.fasta@195.fasta@196.fasta@197.fasta@198.fasta@199.fasta@2.3.fasta@2.fasta@20.fasta@200.fasta@201.fasta@202.fasta@203.fasta@204.fasta@205.fasta@206.fasta@207.fasta@208.fasta@209.fasta@21.fasta@210.fasta@211.fasta@212.fasta@214.fasta@215.fasta@216.fasta@217.fasta@219.fasta@22.fasta@220.fasta@222.fasta@223.fasta@224.fasta@225.fasta@228.fasta@23.fasta@231.fasta@232.fasta@236.fasta@237.fasta@239.fasta@24.fasta@240.fasta@244.fasta@25.fasta@26.fasta@27.fasta@28.fasta@29.fasta@3.fasta@30.fasta@31.fasta@32.fasta@33.fasta@34.fasta@35.fasta@36.fasta@37.fasta@38.fasta@39.fasta@40.fasta@41.fasta@42.fasta@43.fasta@44.fasta@45.fasta@46.fasta@47.fasta@48.fasta@49.fasta@50.fasta@51.fasta@52.fasta@53.fasta@54.fasta@55.fasta@56.fasta@57.fasta@58.fasta@59.fasta@60.fasta@61.fasta@62.fasta@63.fasta@64.fasta@65.fasta@67.fasta@68.fasta@69.fasta@70.fasta@71.fasta@74.fasta@75.fasta@76.fasta@77.fasta@78.fasta@79.fasta@8.fasta@80.fasta@81.fasta@82.fasta@83.fasta@84.fasta@85.fasta@86.fasta@87.fasta@88.fasta@89.fasta@9.fasta@90.fasta@91.fasta@92.fasta@93.fasta@94.fasta@95.fasta@96.fasta@97.fasta@98.fasta@99.fasta@FOC_10.fasta@FOC_3.fasta@FOC_5.1.fasta@FOC_7.1.fasta@NT12001_189.fasta@NT12002_188.fasta@NT12003_214.fasta@NT12004_22.fasta@NT12005_17.fasta@NT12006_215.fasta@NT12007_60.fasta@NT12008_101.fasta@NT12009_154.fasta@NT12010_146.fasta@NT12011_132.fasta@NT12012_88.fasta@NT12013_22.fasta@NT12014_22.fasta@NT12015_11.fasta@NT12016_152.fasta@NT12017_35.fasta@NT12018_216.fasta@NT12019_36.fasta@NT12020_217.fasta@NT12021_218.fasta@NT12022_37.fasta@NT12023_38.fasta@NT12024_39.fasta@NT12025_40.fasta@NT12026_41.fasta@NT12027_219.fasta@NT12028_42.fasta@NT12029_43.fasta@NT12030_220.fasta@NT12031_44.fasta@NT12032_221.fasta@NT12033_45.fasta@NT12034_46.fasta@NT12035_47.fasta@NT12036_48.fasta@NT12037_222.fasta@NT12038_49.fasta@NT12039_223.fasta@NT12040_50.fasta@NT12041_224.fasta@NT12042_225.fasta@NT12043_226.fasta@NT12044_51.fasta@NT12045_52.fasta@NT12046_227.fasta@NT12047_228.fasta@NT12048_229.fasta@NT12049_53.fasta@NT12050.fasta@NT12051_54.fasta@NT12052_231.fasta@NT12053_55.fasta@NT12054_232.fasta@NT12055_233.fasta@NT12056_234.fasta@NT12057_235.fasta@NT12058_236.fasta@NT12059_56.fasta@NT12060_237.fasta@NT12061_238.fasta@NT12062_239.fasta@NT12063_57.fasta@NT12064_748.fasta@NT12065_240.fasta@NT12066_241.fasta@NT12067_242.fasta@NT12068_243.fasta@NT12069_244.fasta@NT12070_245.fasta@NT12071_246.fasta@NT12072_58.fasta@NT12073_247.fasta@NT12074_248.fasta@NT12075_249.fasta@NT12076_250.fasta@NT12077_251.fasta@NT12078.fasta@NT12079_253.fasta@NT12080_254.fasta@NT12081_59.fasta@NT12082_255.fasta@NT12083_256.fasta@NT12084_257.fasta@NT12085_258.fasta@NT12086_259.fasta@NT12087_260.fasta@NT12088_261.fasta@NT12089_32.fasta@NT12090_204.fasta@NT12091_12.fasta@NT12092_10.fasta@NT12093_170.fasta@NT12094_30.fasta@NT12095_187.fasta@NT12096_3.fasta@NT12097_202.fasta@NT12098_101.fasta@NT12099_89.fasta@NT12100_90.fasta@NT12101_91.fasta@NT12102_92.fasta@NT12103_262.fasta@NT12104_93.fasta@NT12105_263.fasta@NT12106_94.fasta@NT12107_264.fasta@NT12108_95.fasta@NT12109_96.fasta@NT12110_97.fasta@NT12111_265.fasta@NT12112_266.fasta@NT12113_98.fasta@NT12114_99.fasta@NT12115_267.fasta@NT12116_100.fasta@NT12117_268.fasta@NT12118_269.fasta@NT12119_102.fasta@NT12120_270.fasta@NT12121_103.fasta@NT12122_271.fasta@NT12123_272.fasta@NT12124_273.fasta@NT12125_274.fasta@NT12126_104.fasta@NT12127_105.fasta@NT12128_275.fasta@NT12129.fasta@NT12130_106.fasta@NT12131_107.fasta@NT12132_277.fasta@NT12133_278.fasta@NT12134_279.fasta@NT12135_280.fasta@NT12136_108.fasta@NT12137_109.fasta@NT12138_281.fasta@NT12140_283.fasta@NT12141_284.fasta@NT12142.fasta@NT12143_286.fasta@NT12144_110.fasta@NT12145_111.fasta@NT12146_287.fasta@NT12147_288.fasta@NT12148_289.fasta@NT12149.fasta@NT12150_291.fasta@NT12151_112.fasta@NT12152_292.fasta@NT12153_293.fasta@NT12154_113.fasta@NT12155_294.fasta@NT12156_295.fasta@NT12157_114.fasta@NT12158_296.fasta@NT12159_750.fasta@NT12160_115.fasta@NT12161_297.fasta@NT12162_298.fasta@NT12163_299.fasta@NT12164_300.fasta@NT12165_301.fasta@NT12166_302.fasta@NT12167_303.fasta@NT12168_304.fasta@NT12169_305.fasta@NT12170_306.fasta@NT12171_307.fasta@NT12172_308.fasta@NT12173_309.fasta@NT12174_310.fasta@NT12175_311.fasta@NT12176_312.fasta@NT12177_313.fasta@NT12178_314.fasta@NT12179_315.fasta@NT12191_751.fasta@NT12192_27.fasta@NT12193_86.fasta@NT12194_61.fasta@NT12195_148.fasta@NT12196_324.fasta@NT12197_151.fasta@NT12198.fasta@NT12199_325.fasta@NT12200.fasta@NT12201.fasta@NT12202_326.fasta@NT12203.fasta@NT12204_755.fasta@NT12205_756.fasta@NT12206.fasta@NT12207.fasta@NT12208.fasta@NT12209.fasta@NT12210_329.fasta@NT12211_330.fasta@NT12212_331.fasta@NT12213_332.fasta@NT12214.fasta@NT12215_333.fasta@NT12216.fasta@NT12217_334.fasta@NT12218_335.fasta@NT12219_757.fasta@NT12220_758.fasta@NT12221_759.fasta@NT12222.fasta@NT12223_336.fasta@NT12224_760.fasta@NT12225.fasta@NT12226_337.fasta@NT12227.fasta@NT12228.fasta@NT12229.fasta@NT12230_338.fasta@NT12231.fasta@NT12232_339.fasta@NT12233_340.fasta@NT12234.fasta@NT12235.fasta@NT12236.fasta@NT12237.fasta@NT12238.fasta@NT12239.fasta@NT12240.fasta@NT12241.fasta@NT12242.fasta@NT12243.fasta@NT12244.fasta@NT12245.fasta@NT12246.fasta@NT12247.fasta@NT12248.fasta@NT12249.fasta@NT12250.fasta@NT12251.fasta@NT12252.fasta@NT12253.fasta@NT12254.fasta@NT12255.fasta@NT12256.fasta@NT12257_761.fasta@NT12258.fasta@NT12259.fasta@NT12260.fasta@NT12261.fasta@NT12262_762.fasta@NT12263.fasta@NT12264.fasta@NT12265.fasta@NT12266.fasta@NT12267.fasta@NT12268.fasta@NT12269.fasta@NT12270.fasta@NT12271.fasta@NT12272.fasta@NT12273.fasta@NT12274.fasta@NT12275.fasta@NT12276.fasta@NT12277.fasta@NT12278.fasta@NT12279.fasta@NT12280.fasta@NT12281.fasta@NT12282.fasta@NT12283.fasta@NT12284.fasta@NT12285_763.fasta@NT12286_764.fasta@NT12287.fasta@NT12288_765.fasta@NT12289.fasta@NT12290.fasta@NT12291.fasta@NT12292_766.fasta@NT12293.fasta@NT12294_341.fasta@NT12295.fasta@NT12296.fasta@NT12297.fasta@NT12298_342.fasta@NT12299_343.fasta@NT12300.fasta@NT12301.fasta@NT12302_344.fasta@NT12303_345.fasta@NT12304.fasta@NT12305_346.fasta@NT12306_347.fasta@NT12307.fasta@NT12308_348.fasta@NT12309_349.fasta@NT12310_767.fasta@NT12311_350.fasta@NT12312.fasta@NT12313.fasta@NT12314_351.fasta@NT12315_352.fasta@NT12316.fasta@NT12317_353.fasta@NT12318.fasta@NT12319_354.fasta@NT12320_355.fasta@NT12321_356.fasta@NT12322.fasta@NT12323.fasta@NT12324.fasta@NT12325_133.fasta@NT12326_143.fasta@NT12327_134.fasta@NT12328_138.fasta@NT12329_142.fasta@NT12330_144.fasta@NT12331_141.fasta@NT12332_139.fasta@NT12333_140.fasta@NT12334_145.fasta@NT12335_135.fasta@NT12336_136.fasta@NT12337_137.fasta@NT12338_2.fasta@NT12339.fasta@NT12340_27.fasta@NT12341_194.fasta@NT12343_357.fasta@NT12344.fasta@NT12345.fasta@NT12346_358.fasta@NT12347_359.fasta@NT12348.fasta@NT12349.fasta@NT12350.fasta@NT12351.fasta@NT12352.fasta@NT12353_360.fasta@NT12354.fasta@NT12355_361.fasta@NT12356.fasta@NT12357.fasta@NT12358_362.fasta@NT12359_363.fasta@NT12360_768.fasta@NT12361_364.fasta@NT12362_365.fasta@NT12363_366.fasta@NT12364_367.fasta@NT12365_368.fasta@NT12367_370.fasta@NT12368.fasta@NT12370_371.fasta@NT12371.fasta@NT12372_372.fasta@NT12373.fasta@NT12375.fasta@NT12376.fasta@NT12377.fasta@NT12378_769.fasta@NT12379.fasta@NT12380.fasta@NT12381.fasta@NT12382_374.fasta@NT12384_770.fasta@NT12385.fasta@NT12386.fasta@NT12387_771.fasta@NT12388_376.fasta@NT12389_772.fasta@NT12390_377.fasta@NT12391_773.fasta@NT12392_378.fasta@NT12393_774.fasta@NT12394_379.fasta@NT12395.fasta@NT12396_380.fasta@NT12397.fasta@NT12398_775.fasta@NT12399.fasta@NT12400_381.fasta@NT12401.fasta@NT12402.fasta@NT12404_383.fasta@NT12405.fasta@NT12406_384.fasta@NT12407_776.fasta@NT12409.fasta@NT12410.fasta@NT12411.fasta@NT12412_386.fasta@NT12413_387.fasta@NT12414_388.fasta@NT12415.fasta@NT12416.fasta@NT12417.fasta@NT12418_389.fasta@NT12419.fasta@NT12420_390.fasta@NT12421_391.fasta@NT12422_392.fasta@NT12423_393.fasta@NT12424_394.fasta@NT12425_395.fasta@NT12426_777.fasta@NT12427_396.fasta@NT12428.fasta@NT12429_397.fasta@NT12430.fasta@NT12432_398.fasta@NT12433.fasta@NT12435_399.fasta@NT12437.fasta@NT12439.fasta@NT12440_401.fasta@NT12441.fasta@NT12442_402.fasta@NT12443.fasta@NT12446.fasta@NT12447_15.fasta@NT12448_16.fasta@NT12449_33.fasta@NT12450_34.fasta@NT12451_76.fasta@NT12452_77.fasta@NT12453_78.fasta@NT12454_79.fasta@NT12457_191.fasta@NT12458_192.fasta@NT12459_193.fasta@NT12460_195.fasta@NT12462_213.fasta@NT12531_172.fasta@NT12532_179.fasta@NT12533_5.fasta@NT12534_659.fasta@NT12535_660.fasta@NT12536_661.fasta@NT12537_662.fasta@NT12538_663.fasta@NT12539_664.fasta@NT12540_665.fasta@NT12542_667.fasta@NT12543_668.fasta@NT12544_669.fasta@NT12545_670.fasta@NT12546_671.fasta@NT12547_672.fasta@NT12548_673.fasta@NT12549_674.fasta@NT12550_675.fasta@NT12551_676.fasta@NT12552_677.fasta@NT12553_678.fasta@NT12554_679.fasta@NT12555_680.fasta@NT12556_681.fasta@NT12557_682.fasta@NT12558_683.fasta@NT12559_684.fasta@NT12560_685.fasta@NT12561_686.fasta@NT12562_687.fasta@NT12563_688.fasta@NT12565_690.fasta@NT12566_691.fasta@NT12567_692.fasta@NT12568_693.fasta@NT12569_694.fasta@NT12570_695.fasta@NT12571_696.fasta@NT12572_697.fasta@NT12573_698.fasta@NT12574_699.fasta@NT12575_700.fasta@NT12576_701.fasta@NT12577_702.fasta@NT12578_703.fasta@NT12579_704.fasta@NT12580_705.fasta@NT12581_706.fasta@NT12582_707.fasta@NT12583_708.fasta@NT12584_709.fasta@NT12585_710.fasta@NT12586_711.fasta@NT12587_712.fasta@NT12588_713.fasta@NT12589_714.fasta@NT12590_715.fasta@NT12591_716.fasta@NT12592_717.fasta@NT12593_718.fasta@NT12594_719.fasta@NT12595_720.fasta@NT12596_721.fasta@NT12597_722.fasta@NT12598_723.fasta@NT12599_724.fasta@NT12600_725.fasta@NT12601_726.fasta@NT12602_727.fasta@NT12603_728.fasta@NT12604_729.fasta@NT12605_730.fasta@NT12606_731.fasta@NT12607_732.fasta@NT12608_733.fasta@NT12609_734.fasta@NT12610_735.fasta@NT12611_736.fasta@NT12612_737.fasta@NT12613.fasta@NT12614_739.fasta@NT12615_740.fasta@NT12616_171.fasta@NT12617_147.fasta@NT12618_181.fasta@NT12619_170.fasta@NT12704_20.fasta@NT12705_21.fasta@Nissle_Tn_parental.fasta@OP50.fasta@SPC_3.1.fasta@SPC_4.fasta


