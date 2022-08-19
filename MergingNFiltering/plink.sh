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
