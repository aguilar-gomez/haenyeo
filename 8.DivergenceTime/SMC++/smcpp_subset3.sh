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
KOR:KOREA1K-313,KOREA1K-83,KOREA1K-36
done
#####################################################################################
nohup bash Kor_vcf2smcS3.sh > outKORconversionSub3 &

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
JEJ:24_S56_L001_chr.bam,25_S57_L001_chr.bam,26_S58_L001_chr.bam
done
#####################################################################################
nohup bash Jej_vcf2smcS3.sh > outJEJconversionSub3 &

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
KOR:KOREA1K-313,KOREA1K-83,KOREA1K-36 \
JEJ:24_S56_L001_chr.bam,25_S57_L001_chr.bam,26_S58_L001_chr.bam
done
#####################################################################################
nohup bash KorJej_vcf2smcS3.sh > outKORJEJconversionSub3 &


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
JEJ:24_S56_L001_chr.bam,25_S57_L001_chr.bam,26_S58_L001_chr.bam \
KOR:KOREA1K-313,KOREA1K-83,KOREA1K-36
done
#####################################################################################
nohup bash JejKor_vcf2smcS3.sh > outJEJKORconversionSub3 &


#Run split to refine the marginal estimates into an estimate of the joint demography:

smc++ split -o split/ KOR/model.final.json JEJ/model.final.json *.smc.gz
smc++ plot allChr_allSamples_jointS3.pdf split/model.final.json
