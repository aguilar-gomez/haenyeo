#PBS scan
/space/s2/diana/korea/haenyeo/Imputation2023/AnalysisNew/outgroupCHS_allSamples/PBS
cp ../allPopsUpdated.fam .

grep "HAE\|JEU\|SEO\|CHB\|KOA" allPopsUpdated.fam |cut -f1|sort|uniq -c
     57 CHB
     27 HAE
     24 JEU
     77 KOA
     26 SEO

grep "HAE\|JEU\|SEO\|CHB\|KOA" allPopsUpdated.fam > samples4scans
cut -f2 samples4scans > inds2keep
sed 's/ /\t/g' pcapops.allchr.maf2.v11.fam > tab.fam


filterbycol.py tab.fam 1 inds2keep inds2keep.fam

Generate vcf
#Extract only Korean and Northern Chinese, keep sites variable in Korean 
plink --file haenyeo.forohanav11 --recode vcf-iid --out selection_scansv11 --keep inds2keep.fam –extract onlyKor.bim
nohup sh extract_for_scan.sh > out.extract4scan &
6377716 variants and 206 people pass filters and QC.


selection_scansv11

At least 75% of the data (this did not work)
vcftools --recode --vcf selection_scansv11.vcf --max-missing .25 --out selection_scansv11_reduced.vcf
After filtering, kept 6377689 out of a possible 6377716 Sites
grep "chr7_155206150" selection_scansv11_reduced.vcf.recode.vcf

Max missing genotypes set to 30 (this is version12)
vcftools --recode --vcf selection_scansv11.vcf --max-missing-count 30 --out selection_scansv11_missingCount
Parameters as interpreted:
    --vcf selection_scansv11.vcf
    --max-missing-count 30
    --out selection_scansv11_missingCount
    --recode
After filtering, kept 5,006,988 out of a possible 6,377,716 Sites

grep "HAE"  samples4scans |cut -f2 >hae.list
grep "JEJ"  samples4scans |cut -f2 >jej.list
grep "SEO\|KOR"  samples4scans |cut -f2 >seo.list
grep "CHB"  samples4scans |cut -f2 >chb.list
grep "HAE\|JEJ"  samples4scans |cut -f2 >jisland.list
cat hae.list jeju.list >haejeju.list

VCFtools (0.1.17)
Warning: Expected at least 2 parts in INFO entry: ID=PR,Number=0,Type=Flag,Description="Provisional reference allele, may not be based on real reference genome">

#A heanyeo+jeju, seul 
vcftools --vcf selection_scansv11.vcf  \
--weir-fst-pop haejeju.list --weir-fst-pop seo.list \
--fst-window-size 20000 --out ./haejeju_seul_20k > fst.w.outA 


#B heanyeo+jeju, chinese
vcftools --vcf selection_scansv11.vcf  \
--weir-fst-pop haejeju.list --weir-fst-pop chb.list \
--fst-window-size 20000 --out ./haejeju_chinese_20k > fst.w.outB 

#C seul, chinese
vcftools --vcf selection_scansv11.vcf  \
--weir-fst-pop seo.list --weir-fst-pop chb.list \
--fst-window-size 20000 --out ./seoul_chinese_20k > fst.w.outC 

By site
#A heanyeo+jeju, seul 
vcftools --vcf selection_scansv11.vcf  \
--weir-fst-pop haejeju.list --weir-fst-pop seo.list \
--out ./haejeju_seoul > fst.outA 

#B heanyeo+jeju, chinese
vcftools --vcf selection_scansv11.vcf  \
--weir-fst-pop haejeju.list --weir-fst-pop chb.list \
--out ./haejeju_chinese > fst.outB 

#C seul, chinese
vcftools --vcf selection_scansv11.vcf  \
--weir-fst-pop seo.list --weir-fst-pop chb.list \
--out ./seoul_chinese > fst.outC 

paste -d "\t" haejeju_seoul.weir.fst  haejeju_chinese.weir.fst seoul_chinese.weir.fst | cut -f1-3,6,9 |grep -v "nan" > hsc4pbs

grep -v "nan" hsc.fst > hsc4pbs


bedtools intersect -wa -wb -a HAE_JEJU_pbs_outliers_top_persitepbs1_25.bed -b annotation_grch38.bed |cut -f1-4,8>HAE_JEJU_pbs_anno25.bed


bedtools intersect -wa -wb -a HAE_JEJU_pbs_outliers_top_persite_chr6pbs1_25.bed -b annotation_grch38.bed |cut -f1-4,8>HAE_JEJU_pbs_anno25_chr6.bed






#By site
#A heanyeo+jeju, seul
vcftools --vcf selection_scansv11_missingCount.recode.vcf \
--weir-fst-pop haejeju.list --weir-fst-pop seo.list \
--out ./haejeju_seoulMC > fst.outAMC

#B heanyeo+jeju, chinese
vcftools --vcf selection_scansv11_missingCount.recode.vcf \
--weir-fst-pop haejeju.list --weir-fst-pop chb.list \
--out ./haejeju_chineseMC > fst.outBMC

#C seul, chinese
vcftools --vcf selection_scansv11_missingCount.recode.vcf \
--weir-fst-pop seo.list --weir-fst-pop chb.list \
--out ./seoul_chineseMC > fst.outCMC

After filtering, kept 5006988 out of a possible 5006988 Sites

nohup sh persite_fst_missingCount.sh >out.missing.fst &
paste -d "\t"  *MC*fst | cut -f1-3,6,9 |grep -v "nan" > hsc4pbsMC

paste -d "\t"  haejeju_seoulMC.weir.fst haejeju_chineseMC.weir.fst  seoul_chineseMC.weir.fst | cut -f1-3,6,9 |grep -v "nan" > hsc4pbsMC

haejeju_seoulMC.weir.fst haejeju_chineseMC.weir.fst  seoul_chineseMC.weir.fst
bedtools intersect -wa -wb -a HAE_JEJU_pbs_outliers_top_persitepbs1_23.bed -b annotation_grch38.bed |cut -f1-4,8>HAE_JEJU_pbs_anno23.bed

grep "chr6_26746360" selection_scansv11_missingCount.recode.vcf 

HAE_JEJU_pbs_outliers_top_persitepbs1_3.48.bed

bedtools intersect -wa -wb -a HAE_JEJU-SEOULfst_outliers_top_persitefst01_25.bed-b annotation_grch38.bed |cut -f1-4,8>HAE_JEJU-SEOUL_fst_anno25.bed






bedtools intersect -wa -wb -a HAE_JEJU_pbs_outliers_top_persitepbs2_2.bed -b annotation_grch38.bed |cut -f1-4,8>HAE_JEJUpbs2_anno2.bed
