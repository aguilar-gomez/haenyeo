#####################################################################################
#Kor_vcf2smc.sh
#Convert vcf for each chromosome
D1=KOREA1K-313
D2=KOREA1K-83
for i in {1..22}; 
do
chr=chr$i
echo $chr
smc++ vcf2smc -d ${D1} ${D2} KoreanJeju.vcf.gz KOR.${D1}.${D2}.$chr.smc.gz $chr \
KOR:KOREA1K-313,KOREA1K-83,KOREA1K-36,KOREA1K-224,62_S98_L002_chr.bam,65_S100_L002_chr.bam,71_S105_L002_chr.bam,76_S108_L002_chr.bam
done
#####################################################################################
nohup bash Kor_vcf2smcS2.sh > outKORconversionSub2 &

smc++ estimate -o KOR 1.25e-8 -p 0.5 KOR.*.smc.gz

#####################################################################################
#Jej_vcf2smc.sh
#Convert vcf for each chromosome
D1=24_S56_L001_chr.bam
D2=25_S57_L001_chr.bam
for i in {1..22}; 
do
chr=chr$i
echo $chr
smc++ vcf2smc -d ${D1} ${D2} KoreanJeju.vcf.gz JEJ.${D1}.${D2}.$chr.smc.gz $chr \
JEJ:24_S56_L001_chr.bam,25_S57_L001_chr.bam,26_S58_L001_chr.bam,27_S59_L001_chr.bam,30_S61_L001_chr.bam,31_S62_L001_chr.bam,41_S72_L001_chr.bam,42_S73_L001_chr.bam
done
#####################################################################################
nohup bash Jej_vcf2smcS2.sh > outJEJconversionSub2 &

smc++ estimate -o JEJ 1.25e-8 -p 0.5 JEJ.*.smc.gz

#Create datasets containing the joint frequency spectrum for both populations:
#####################################################################################
#KorJej_vcf2smc.sh
#Convert vcf for each chromosome
D1=KOREA1K-313
D2=24_S56_L001_chr.bam
for i in {1..22}; 
do
chr=chr$i
echo $chr
smc++ vcf2smc -d ${D1} ${D2} KoreanJeju.vcf.gz KORJEJ.${D1}.${D2}.$chr.smc.gz $chr \
KOR:KOREA1K-313,KOREA1K-83,KOREA1K-36,KOREA1K-224,62_S98_L002_chr.bam,65_S100_L002_chr.bam,71_S105_L002_chr.bam,76_S108_L002_chr.bam \
JEJ:24_S56_L001_chr.bam,25_S57_L001_chr.bam,26_S58_L001_chr.bam,27_S59_L001_chr.bam,30_S61_L001_chr.bam,31_S62_L001_chr.bam,41_S72_L001_chr.bam,42_S73_L001_chr.bam
done
#####################################################################################
nohup bash KorJej_vcf2smcS2.sh > outKORJEJconversionSub &


#####################################################################################
#JejKor_vcf2smc.sh
#Convert vcf for each chromosome
D1=24_S56_L001_chr.bam
D2=KOREA1K-313
for i in {1..22}; 
do
chr=chr$i
echo $chr
smc++ vcf2smc -d ${D1} ${D2} KoreanJeju.vcf.gz JEJKOR.${D1}.${D2}.$chr.smc.gz $chr \
JEJ:24_S56_L001_chr.bam,25_S57_L001_chr.bam,26_S58_L001_chr.bam,27_S59_L001_chr.bam,30_S61_L001_chr.bam,31_S62_L001_chr.bam,41_S72_L001_chr.bam,42_S73_L001_chr.bam \
KOR:KOREA1K-313,KOREA1K-83,KOREA1K-36,KOREA1K-224,62_S98_L002_chr.bam,65_S100_L002_chr.bam,71_S105_L002_chr.bam,76_S108_L002_chr.bam
done
#####################################################################################
nohup bash JejKor_vcf2smcS2.sh > outJEJKORconversionSub &


#Run split to refine the marginal estimates into an estimate of the joint demography:

smc++ split -o split/ KOR/model.final.json JEJ/model.final.json *.smc.gz
smc++ plot allChr_allSamples_joint.pdf split/model.final.json
