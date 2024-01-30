#####################################################################################
#Kor_vcf2smc.sh
#Convert vcf for each chromosome
D1=KOREA1K-366
D2=KOREA1K-175
for i in {1..22}; 
do
chr=chr$i
echo $chr
smc++ vcf2smc -d ${D1} ${D2} KoreanJeju.vcf.gz KOR.${D1}.${D2}.$chr.smc.gz $chr \
KOR:KOREA1K-175,KOREA1K-366,KOREA1K-45,KOREA1K-70,KOREA1K-97,KOREA1K-30,KOREA1K-217,\
KOREA1K-7,KOREA1K-231,KOREA1K-92,KOREA1K-360,KOREA1K-50,KOREA1K-167,KOREA1K-63,KOREA1K-105,\
KOREA1K-126,KOREA1K-84,KOREA1K-106,KOREA1K-268,KOREA1K-325,KOREA1K-144,KOREA1K-87,KOREA1K-211,\
KOREA1K-302,KOREA1K-313,KOREA1K-83,KOREA1K-36,KOREA1K-224,62_S98_L002_chr.bam,65_S100_L002_chr.bam,\
71_S105_L002_chr.bam,76_S108_L002_chr.bam,79_S111_L002_chr.bam,89_S177_L003_chr.bam,\
90_S178_L003_chr.bam,KOREA1K-194,83_S173_L003_chr.bam
done
#####################################################################################
nohup bash Kor_vcf2smc.sh > outKORconversionUA &

smc++ estimate -o KOR 1.25e-8 -p 0.5 --thinning 400 --cores 8 KOR.*.smc.gz

#####################################################################################
#Jej_vcf2smc.sh
#Convert vcf for each chromosome
D1=10_S44_L001_chr.bam
D2=1_S35_L001_chr.bam
for i in {1..22}; 
do
chr=chr$i
echo $chr
smc++ vcf2smc -d ${D1} ${D2} KoreanJeju.vcf.gz JEJ.${D1}.${D2}.$chr.smc.gz $chr \
JEJ:KOREA1K-248,10_S44_L001_chr.bam,13_S47_L001_chr.bam,14_S48_L001_chr.bam,15_S49_L001_chr.bam,\
18_S51_L001_chr.bam,19_S52_L001_chr.bam,1_S35_L001_chr.bam,20_S53_L001_chr.bam,23_S55_L001_chr.bam,\
24_S56_L001_chr.bam,25_S57_L001_chr.bam,26_S58_L001_chr.bam,27_S59_L001_chr.bam,30_S61_L001_chr.bam,\
31_S62_L001_chr.bam,41_S72_L001_chr.bam,42_S73_L001_chr.bam,45_S76_L002_chr.bam,49_S80_L002_chr.bam,\
4_S38_L001_chr.bam,51_S82_L002_chr.bam,5_S39_L001_chr.bam,6_S40_L001_chr.bam,7_S41_L001_chr.bam,\
8_S42_L001_chr.bam,DNA_32_S63_L001_chr.bam
done
#####################################################################################
nohup bash Jej_vcf2smc.sh > outJEJconversionUA &

smc++ estimate -o JEJ 1.25e-8 -p 0.5 --thinning 400 --cores 8 JEJ.*.smc.gz

#Create datasets containing the joint frequency spectrum for both populations:
#####################################################################################
#KorJej_vcf2smc.sh
#Convert vcf for each chromosome
D1=KOREA1K-366
D2=10_S44_L001_chr.bam
for i in {1..22}; 
do
chr=chr$i
echo $chr
smc++ vcf2smc -d ${D1} ${D2} KoreanJeju.vcf.gz KORJEJ.${D1}.${D2}.$chr.smc.gz $chr \
KOR:KOREA1K-175,KOREA1K-366,KOREA1K-45,KOREA1K-70,KOREA1K-97,KOREA1K-30,KOREA1K-217,\
KOREA1K-7,KOREA1K-231,KOREA1K-92,KOREA1K-360,KOREA1K-50,KOREA1K-167,KOREA1K-63,KOREA1K-105,\
KOREA1K-126,KOREA1K-84,KOREA1K-106,KOREA1K-268,KOREA1K-325,KOREA1K-144,KOREA1K-87,KOREA1K-211,\
KOREA1K-302,KOREA1K-313,KOREA1K-83,KOREA1K-36,KOREA1K-224,62_S98_L002_chr.bam,65_S100_L002_chr.bam,\
71_S105_L002_chr.bam,76_S108_L002_chr.bam,79_S111_L002_chr.bam,89_S177_L003_chr.bam,\
90_S178_L003_chr.bam,KOREA1K-194,83_S173_L003_chr.bam \
JEJ:KOREA1K-248,10_S44_L001_chr.bam,13_S47_L001_chr.bam,14_S48_L001_chr.bam,15_S49_L001_chr.bam,\
18_S51_L001_chr.bam,19_S52_L001_chr.bam,1_S35_L001_chr.bam,20_S53_L001_chr.bam,23_S55_L001_chr.bam,\
24_S56_L001_chr.bam,25_S57_L001_chr.bam,26_S58_L001_chr.bam,27_S59_L001_chr.bam,30_S61_L001_chr.bam,\
31_S62_L001_chr.bam,41_S72_L001_chr.bam,42_S73_L001_chr.bam,45_S76_L002_chr.bam,49_S80_L002_chr.bam,\
4_S38_L001_chr.bam,51_S82_L002_chr.bam,5_S39_L001_chr.bam,6_S40_L001_chr.bam,7_S41_L001_chr.bam,\
8_S42_L001_chr.bam,DNA_32_S63_L001_chr.bam
done
#####################################################################################
nohup bash KorJej_vcf2smc.sh > outKORJEJconversionUA &


#####################################################################################
#JejKor_vcf2smc.sh
#Convert vcf for each chromosome
D1=10_S44_L001_chr.bam
D2=KOREA1K-366
for i in {1..22}; 
do
chr=chr$i
echo $chr
smc++ vcf2smc -d ${D1} ${D2} KoreanJeju.vcf.gz JEJKOR.${D1}.${D2}.$chr.smc.gz $chr \
JEJ:KOREA1K-248,10_S44_L001_chr.bam,13_S47_L001_chr.bam,14_S48_L001_chr.bam,15_S49_L001_chr.bam,\
18_S51_L001_chr.bam,19_S52_L001_chr.bam,1_S35_L001_chr.bam,20_S53_L001_chr.bam,23_S55_L001_chr.bam,\
24_S56_L001_chr.bam,25_S57_L001_chr.bam,26_S58_L001_chr.bam,27_S59_L001_chr.bam,30_S61_L001_chr.bam,\
31_S62_L001_chr.bam,41_S72_L001_chr.bam,42_S73_L001_chr.bam,45_S76_L002_chr.bam,49_S80_L002_chr.bam,\
4_S38_L001_chr.bam,51_S82_L002_chr.bam,5_S39_L001_chr.bam,6_S40_L001_chr.bam,7_S41_L001_chr.bam,\
8_S42_L001_chr.bam,DNA_32_S63_L001_chr.bam \
KOR:KOREA1K-175,KOREA1K-366,KOREA1K-45,KOREA1K-70,KOREA1K-97,KOREA1K-30,KOREA1K-217,\
KOREA1K-7,KOREA1K-231,KOREA1K-92,KOREA1K-360,KOREA1K-50,KOREA1K-167,KOREA1K-63,KOREA1K-105,\
KOREA1K-126,KOREA1K-84,KOREA1K-106,KOREA1K-268,KOREA1K-325,KOREA1K-144,KOREA1K-87,KOREA1K-211,\
KOREA1K-302,KOREA1K-313,KOREA1K-83,KOREA1K-36,KOREA1K-224,62_S98_L002_chr.bam,65_S100_L002_chr.bam,\
71_S105_L002_chr.bam,76_S108_L002_chr.bam,79_S111_L002_chr.bam,89_S177_L003_chr.bam,\
90_S178_L003_chr.bam,KOREA1K-194,83_S173_L003_chr.bam 
done
#####################################################################################
nohup bash JejKor_vcf2smc.sh > outJEJKORconversionUA &


#Run split to refine the marginal estimates into an estimate of the joint demography:

smc++ split --thinning 800 --cores 8 -o split/ KOR/model.final.json JEJ/model.final.json *.smc.gz
smc++ plot allChr_Unrelated_notAdmixed_jointEstimation.pdf split/model.final.json
