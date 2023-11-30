#Download reference panel from 1000 genomes
#Reference panel
#!/bin/bash
for i in {1..22}; 
do echo $i 
chr=chr$i
echo $chr
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.${chr}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz{,.tbi};
done


#Only keep female individuals. From all populations, I am ONLY using female individuals.
#Download list of samples from 100 genes and grep "female" and populations of interest
cat *Phase3|grep "female" > AsianFemale_Phase3_list
cut -f4 AsianFemale_Phase3_list |uniq -c

#Extract names of samples for bcftools
cut -f1  Female_Phase3_list> indvsamples1000genomes 

for i in {1..22};
do echo $i
chr=chr$i
echo $chr
bcftools view -S indvsamples1000genomes -o $chr.female1000genomes.vcf.gz --force-samples -Oz ALL.${chr}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz &
done


#Get Strict genome accessibility mask
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20160622_genome_mask_GRCh38/StrictMask/20160622.allChr.mask.bed

#Merge vcfs and add mask

for chr in chr{1..22}
do  
echo $chr
bcftools index ${chr}.allK2023Mask.vcf.gz
bcftools index $chr.female1000genomes.vcf.gz
bcftools merge ${chr}.allK2023Mask.vcf.gz $chr.female1000genomes.vcf.gz --threads 10 -O z -o ${chr}.allSamples2023.vcf.gz
bcftools index ${chr}.allSamples2023.vcf.gz
bcftools view -m2 -M2 -v snps -R 20160622.allChr.mask.bed -o  ${chr}.allSamplesMask2023.vcf.gz  -O z ${chr}.allSamples2023.vcf.gz –threads 15 -S ^samples2exclude
done

#Apply extra filters

for chr in chr{1..22}
do  
echo $chr
bcftools +fill-tags ${chr}.allSamplesMask2023.vcf.gz -- -t all |grep -v 'INFO=REF_SWITCH'|bcftools view -e'ExcHet<=.000001 || MAF=0' -Oz  --threads 5 -o $chr.filteredAlll.vcf.gz
done

bcftools query -l $chr.filteredAlll.vcf.gz


#Maf filter and no missing data
for chr in chr{1..22}
do  
echo $chr
bcftools view -e 'GT[*] = "mis"' -q 0.01:minor $chr.filteredAlll.vcf.gz -Oz -o $chr.filteredAllMaf1.vcf.gz
done


