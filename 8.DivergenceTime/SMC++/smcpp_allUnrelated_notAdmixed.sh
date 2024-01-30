#####################################################################################
#Kor_vcf2smc.sh
#Convert vcf for each chromosome
D1=KOREA1K-348
D2=KOREA1K-175
for i in {1..22}; 
do
chr=chr$i
echo $chr
smc++ vcf2smc -d ${D1} ${D2} KoreanJeju.vcf.gz KOR.${D1}.${D2}.$chr.smc.gz $chr \
KOR:KOREA1K-348,KOREA1K-175,KOREA1K-366,KOREA1K-45,KOREA1K-70,KOREA1K-97,KOREA1K-30,KOREA1K-217
done
#####################################################################################
nohup bash Kor_vcf2smc.sh > outKORconversionSub &

smc++ estimate -o KOR 1.25e-8 -p 0.5 KOR.*.smc.gz

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
JEJ:10_S44_L001_chr.bam,1_S35_L001_chr.bam,20_S53_L001_chr.bam
done
#####################################################################################
nohup bash Jej_vcf2smc.sh > outJEJconversionSub &

smc++ estimate -o JEJ 1.25e-8 -p 0.5 JEJ.*.smc.gz

#Create datasets containing the joint frequency spectrum for both populations:
#####################################################################################
#KorJej_vcf2smc.sh
#Convert vcf for each chromosome
D1=KOREA1K-348
D2=10_S44_L001_chr.bam
for i in {1..22}; 
do
chr=chr$i
echo $chr
smc++ vcf2smc -d ${D1} ${D2} KoreanJeju.vcf.gz KORJEJ.${D1}.${D2}.$chr.smc.gz $chr \
KOR:KOREA1K-348,KOREA1K-175,KOREA1K-366,KOREA1K-45,KOREA1K-70,KOREA1K-97,KOREA1K-30,KOREA1K-217 \
JEJ:10_S44_L001_chr.bam,1_S35_L001_chr.bam,20_S53_L001_chr.bam
done
#####################################################################################
nohup bash KorJej_vcf2smc.sh > outKORJEJconversionSub &


#####################################################################################
#JejKor_vcf2smc.sh
#Convert vcf for each chromosome
D1=10_S44_L001_chr.bam
D2=KOREA1K-348
for i in {1..22}; 
do
chr=chr$i
echo $chr
smc++ vcf2smc -d ${D1} ${D2} KoreanJeju.vcf.gz JEJKOR.${D1}.${D2}.$chr.smc.gz $chr \
JEJ:10_S44_L001_chr.bam,1_S35_L001_chr.bam,20_S53_L001_chr.bam \
KOR:KOREA1K-348,KOREA1K-175,KOREA1K-366,KOREA1K-45,KOREA1K-70,KOREA1K-97,KOREA1K-30,KOREA1K-217 

done
#####################################################################################
nohup bash JejKor_vcf2smc.sh > outJEJKORconversionSub &


#Run split to refine the marginal estimates into an estimate of the joint demography:

smc++ split -o split/ KOR/model.final.json JEJ/model.final.json *.smc.gz
smc++ plot allChr_allSamples_joint.pdf split/model.final.json
