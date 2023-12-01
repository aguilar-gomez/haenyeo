#Calculate pi using vcftools diveristy
grep "JEJ" popAssignment |cut -f1 > JejuPop
grep "KOR" popAssignment |cut -f1 > KorPop

bcftools view -S JejuPop KoreanJeju.vcf.gz |grep -v -f Flipped2exclude|gzip > Jeju.vcf.gz
bcftools view -S KorPop KoreanJeju.vcf.gz |grep -v -f Flipped2exclude|gzip > Kor.vcf.gz


#Set window size larger than chr1 size to calculate pi per chromosome
vcftools --gzvcf Jeju.vcf.gz --window-pi 250000000 --out Jeju_pi
awk -v OFS='\t' 'NR==1 {print $0, "POPSIZE"} NR>1 { $6 = $5 / (4 * 1.25 * 10^(-8)); print }' Jeju_pi.windowed.pi > Jeju_estPopSize_fromPi

vcftools --gzvcf Kor.vcf.gz --window-pi 250000000 --out Kor_pi
awk -v OFS='\t' 'NR==1 {print $0, "POPSIZE"} NR>1 { $6 = $5 / (4 * 1.25 * 10^(-8)); print }' Kor_pi.windowed.pi > Kor_estPopSize_fromPi



#Extract labels we had from OHANA
grep "JEJ" popAssignment |cut -f1 > JejuPop
grep "KOR" popAssignment |cut -f1 > KorPop
cat JejuPop KorPop > allKorean

#We used chromosome one to estimate divergence time
chr=1
#The original file has the related files that we removed
#extract KOR and JEJ samples 
bcftools view -S allKorean chr1.allK2023.vcf.gz >chr1.subsetKorean.vcf

#Generate Relate file format
PATH_TO_RELATE=~/programs/relate_ARGS/
$PATH_TO_RELATE/bin/RelateFileFormats \
                --mode ConvertFromVcf  \
                --haps allpos.haps  	\
                --sample allpos.sample 	\
                -i  chr1.subsetKorean 

#Generate the ARGs
#Population size 
#N=2Ne
#Ne calculated using pi (calculated above using vcftools and awk), the average across chromosomes was ~5000
#theta=4Nemu, theta~pi
#Ne=pi/(4mu)
$PATH_TO_RELATE/bin/Relate \
      --mode All \
      -m 1.25e-8 \
      -N 10000 \
      --haps allpos.haps \
      --sample allpos.sample \
      --map mapchr1.txt \
      --seed 1 \
      -o allpops


#Make sure the pop labels are in the correct order
import pandas as pd

popA=pd.read_table("popAssignment",header=None)
popK=pd.read_table("allKorean",header=None)
df=popK.merge(popA)

forRelate=pd.DataFrame({'sample': df[0],
	                      'population': df[1],
	                      'group': df[1],
	                      'sex': 'NA'})

forRelate.to_csv("allSamplesOrder.popfile",sep="\t",index=False)

#Jointly Estimate Population Size of JEJ and KOR
PATH_TO_RELATE=~/programs/relate_ARGS/
$PATH_TO_RELATE/scripts/EstimatePopulationSize/EstimatePopulationSize.sh \
              -i allpops \
              --threads 10 \
              -m 1.25e-8 \
              --years_per_gen 25
              --seed 1\
              --poplabels allSamplesOrder.popfile\
              -o gentime25 \




