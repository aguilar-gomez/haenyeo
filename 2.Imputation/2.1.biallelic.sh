for file in *imputed.vcf.gz
do
chr=`echo $file | cut -d "." -f 1`
echo $chr
bcftools view -G -m 2 -M 2 -v snps $file -Oz -o $chr.biallelic.sites.vcf.gz
bcftools index -f $chr.biallelic.sites.vcf.gz	
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' $chr.biallelic.sites.vcf.gz	 | bgzip -c > $chr.biallelic.sites.tsv.gz	
tabix -s1 -b2 -e2 $chr.biallelic.sites.tsv.gz 
done


