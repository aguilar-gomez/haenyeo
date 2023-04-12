GLIMPSE=~/programs/GLIMPSE-1.1.0/static_bins
for bcf in ../Ligated/*merged.bcf
do
name=$(echo $bcf|cut -d"/" -f3 |cut -d"." -f1-2); 
echo $name
VCF=$bcf
OUT=$name.phased.bcf
$GLIMPSE/GLIMPSE_sample_static --input ${VCF} --solve --output ${OUT}
bcftools index -f ${OUT} 	
done
