#Convert to plink
#PLINK v1.90b6.7
for chr in chr{1..22}
do  
echo $chr
plink --vcf  $chr.hetfil.vcf.gz --make-bed --out snpspopsv11.$chr --threads 10 --const-fid 0 --vcf-half-call m --set-missing-var-ids @:#$1,$2 --snps-only --biallelic-only
done

#Merge all chromosomes
#Remove unwanted individuals (i.e. related individuals)
plink --bfile snpspopsv11.chr1 --merge-list file2merge --make-bed --out snpversion11.allchr --remove samples2removev11

#Principal Component Analysis
plink --bfile snpversion11.allchr --pca --out pcapops.allchr.maf2.v11 --threads 10 --maf 0.02 --make-bed 

#Prepare plink for OHANA
plink --bfile pcapops.allchr.maf2.v11 --recode 12 tab --out haenyeo.forohanav11 

#Extract only Korean and Northern Chinese, keep sites variable in Korean 
plink --file haenyeo.forohanav11 --recode vcf-iid --out selection_scansv11 --keep inds2keep.fam –extract onlyKor.bim

#Remove positions with missing data: max missing genotypes set to 30 in Korean populations, this does not keep all the headers and flags
vcftools --recode --vcf selection_scansv11.vcf --max-missing-count 30 --out selection_scansv11_missingCount

#Merge unfiltered vcfs withh all data
bcftools concat chr1.filter11.vcf.gz chr2.filter11.vcf.gz chr3.filter11.vcf.gz chr4.filter11.vcf.gz chr5.filter11.vcf.gz chr6.filter11.vcf.gz chr7.filter11.vcf.gz chr8.filter11.vcf.gz chr9.filter11.vcf.gz chr10.filter11.vcf.gz chr11.filter11.vcf.gz chr12.filter11.vcf.gz chr13.filter11.vcf.gz chr14.filter11.vcf.gz chr15.filter11.vcf.gz chr16.filter11.vcf.gz chr17.filter11.vcf.gz chr18.filter11.vcf.gz chr19.filter11.vcf.gz chr20.filter11.vcf.gz chr21.filter11.vcf.gz chr22.filter11.vcf.gz -o allchr.version11.vcf.gz -Oz

#Extract positions from file with missing counts filter and format them (add chr)
bcftools query -f 'chr%CHROM\t%POS\n' selection_scansv11_missingCount.recode.vcf.gz > missingCount.pos_version12_chr &

#Generate new version with last filter and preserving all flags
nohup bcftools view -R missingCount.pos_version12_chr allchr.version11.vcf.gz -Oz -o allchr.version12.vcf.gz &

#We decided to remove the positions with a switch in reference allele because we cannot be sure the fixed differences are not due to that and the proportion of SNPs removed is very small. Other positions in LD with them will give us the results if they should be kept
bcftools filter -e 'INFO=REF_SWITCH' allchr.version12.vcf.gz -Oz -o allchr.version13.vcf.gz

#bcftools filter did not remove all, use grep
zcat allchr.version13.vcf.gz | grep -v "REF_SWITCH"> allchr.version13.grep.vcf &

#Convert to plink again 
plink --vcf allchr.version13.grep.vcf --make-bed --out allchr.version14 --remove samples2removev11 --const-fid 0 --vcf-half-call m
 
#Principal Component Analysis
plink --bfile allchr.version14 --pca --threads 10

