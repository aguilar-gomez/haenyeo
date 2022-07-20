#Use this: https://odelaneau.github.io/GLIMPSE/

#rename_bam.sh
#Rename chr in bams:
for file in *38.bam
do 
  filename=`echo $file | cut -d "." -f 1`
  samtools view -H $file | samtools view -H $file | awk -F '\t' '{for(x=1;x<=NF;x++)if($x~/CM.*?/){gsub(/CM.*?L/,"chr"++i"\tL")}}1'  | samtools reheader - $file > ${filename}_chr.bam
  samtools index ${filename}_chr.bam
  rm $file
done

nohup sh rename_bam.sh >out.rename &

#Generate haplotype files, the Reference panel is data we got from the Korea 1K project (400 random individuals from their datset)
#Follow GLIMPSE tutorial
for chr in chr{1..22}
do
bcftools view -Ob -o $chr.imputed.bcf $chr.imputed.vcf.gz --threads 5 
bcftools index -f $chr.imputed.bcf 
done

#biallelic.sh
#Extract only biallelic positions, compress vcf files, generate tsv and indexes
for file in *imputed.vcf.gz
do
chr=`echo $file | cut -d "." -f 1`
echo $chr
bcftools view -G -m 2 -M 2 -v snps $file -Oz -o $chr.biallelic.sites.vcf.gz
bcftools index -f $chr.biallelic.sites.vcf.gz    
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' $chr.biallelic.sites.vcf.gz     | bgzip -c > $chr.biallelic.sites.tsv.gz    
tabix -s1 -b2 -e2 $chr.biallelic.sites.tsv.gz 
done

nohup sh biallelic.sh >out.biallelic &


#Computing GLs for a single individual at specific positions
#As an output of this step we have a VCF file format containing genotype likelihoods at each site in the reference panel.
#call.likelihood.sh
BAM=$1
for CHR in chr{1..22}
do 
echo $CHR
indv=`echo $BAM | cut -d "_" -f 1-2`
echo $indv $CHR
VCF=/media/external3/diana/reference_korean/$CHR.biallelic.sites.vcf.gz
TSV=/media/external3/diana/reference_korean/$CHR.biallelic.sites.tsv.gz
REFGEN=/media/external3/diana/reference/GRCh38/grch38.fasta
OUT=$indv.$CHR.vcf.gz
bcftools mpileup -f ${REFGEN} -I -E -a 'FORMAT/DP' -T ${VCF} -r $CHR ${BAM} -Ou | bcftools call -Aim -C alleles -T ${TSV} -Oz -o ${OUT}
bcftools index -f ${OUT}
done


#Split the genome into chunks
for chr in chr{1..22}
do 
echo $chr
GLIMPSE_chunk_static --input $chr.biallelic.sites.vcf.gz --region $chr --window-size 2000000 --buffer-size 200000 --output chunks.$chr.txt --thread 5
done

#Impute and phase a whole chromosome
VCF=$1  #genotype likelihoods of sample in relevant positions
chr=$2 
name=$(echo $VCF| cut -d"_" -f1-2)
REF=/media/external3/diana/reference_korean/$chr.imputed.bcf 
MAP=/home/diana/programs/GLIMPSE/maps/genetic_maps.b38/$chr.b38.gmap.gz
while IFS="" read -r LINE || [ -n "$LINE" ];
do
ID=$(echo $LINE | cut -d" " -f1)
IRG=$(echo $LINE | cut -d" " -f3)
ORG=$(echo $LINE | cut -d" " -f4)
OUT=Imputed/${name}.${chr}.${ID}.bcf
/home/diana/programs/GLIMPSE/static_bins/GLIMPSE_phase_static --input ${VCF} --reference ${REF} --map ${MAP} --input-region ${IRG} --output-region ${ORG} --output ${OUT}
bcftools index -f ${OUT}
done < /media/external3/diana/reference_korean/chunks.$chr.txt

#Ligate multiple chunks together
VCF=$1  #genotype likelihoods of sample in relevant positions
chr=$2 
name=$(echo $VCF| cut -d"." -f1)
LST=Imputed/$name.list.${chr}.txt
ls Imputed/${name}.${chr}.*.bcf > ${LST}
OUT=Ligated/${name}.${chr}.merged.bcf
/home/diana/programs/GLIMPSE/static_bins/GLIMPSE_ligate_static --input ${LST} --output $OUT
bcftools index -f ${OUT}

#Merge all samples into one file:
for chr in chr{1..22}
do  
echo $chr
bcftools merge -o haenyeo.$chr.vcf.gz -Oz *$chr.*bcf &
done

 


