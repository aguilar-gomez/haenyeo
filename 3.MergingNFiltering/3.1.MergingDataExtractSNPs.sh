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

#Apply mask using option -R (regions)
#Biallelic SNPs
for chr in chr{1..22}
do  
echo $chr
bcftools merge ${chr}.Pops.chr.vcf.gz ${chr}.imputed.vcf.gz --threads 10 -O z -o ${chr}.allsamples.vcf.gz 
bcftools view -m2 -M2 -v snps -R ./masks/20160622.allChr.mask.bed -o $chr.mask.vcf.gz -O z $chr.allsamples.vcf.gz –threads 15
done

#Exclude minor allele frequency to remove samples that are non-variant. 
#Because it is human populations some sites are flagged as SNPs because they are variant in other human populations (not our specific subset)
#Apply excess of heterozygosity filter
for chr in chr{1..22}
do  
echo $chr
bcftools +fill-tags $chr.mask.vcf.gz  -- -t all|bcftools view -e'ExcHet<=.000001'|bcftools view -e'MAF=0' -O z --threads 5 >  $chr.hetfil.vcf.gz
done

