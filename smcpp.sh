

#####################################################################################
#Kor_vcf2smc.sh
#Convert vcf for each chromosome
D1=KOREA1K-348
D2=KOREA1K-175
MASK=20160622.allChr.mask.bed.gz
for i in {1..22}; 
do
chr=chr$i
echo $chr
smc++ vcf2smc -m ${MASK} -d ${D1} ${D2} KoreanJeju.vcf.gz KOR.${D1}.${D2}.$chr.smc.gz $chr \
KOR:KOREA1K-348,KOREA1K-175,KOREA1K-366,KOREA1K-45,KOREA1K-70,KOREA1K-97,KOREA1K-30,KOREA1K-217
done
#####################################################################################
nohup bash Kor_vcf2smc.sh > outKORconversionSub &
