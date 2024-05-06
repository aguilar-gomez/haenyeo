#PBS scan
/space/s2/diana/korea/haenyeo/Imputation2023/AnalysisNew/outgroupCHS_allSamples/PBS
cut -f1 -d " " selscan3popsCHSv2.fam|sort |uniq -c
     53 CHS
     51 JEJ
    108 KOR

#Generate vcf
#Extract only Korean and Northern Chinese, keep sites variable in Korean 
plink --file selscan3popsCHSv2 --recode vcf-iid --out selection_scans 
5673091 variants and 212 people pass filters and QC.

grep "JEJ" selscan3popsCHSv2.fam |cut -f2 -d " " > jej.list
grep "KOR" selscan3popsCHSv2.fam |cut -f2 -d " " > kor.list
grep "CHS" selscan3popsCHSv2.fam |cut -f2 -d " " > chs.list

#A jeju, seoul 
vcftools --vcf selection_scans.vcf  \
--weir-fst-pop jej.list --weir-fst-pop kor.list \
--out ./jej_kor > fst.outA 

#B jeju, chinese
vcftools --vcf selection_scans.vcf  \
--weir-fst-pop jej.list --weir-fst-pop chs.list \
--out ./jej_chinese > fst.outB 

#C seul, chinese
vcftools --vcf selection_scans.vcf  \
--weir-fst-pop kor.list --weir-fst-pop chs.list \
--out ./kor_chinese > fst.outC 

paste -d "\t" jej_kor.weir.fst  jej_chinese.weir.fst kor_chinese.weir.fst | cut -f1-3,6,9 |grep -v "nan" > JEJKORCHB4pbs

#Run PBS script
Rscript PBSfromvcftools.R titleOfPlot pop0 pop1 pop2 outquantile color pop012pbs name

#Delete below
bedtools intersect -wa -wb -a HAE_JEJU_pbs_outliers_top_persitepbs1_25.bed -b annotation_grch38.bed |cut -f1-4,8>HAE_JEJU_pbs_anno25.bed
bedtools intersect -wa -wb -a HAE_JEJU_pbs_outliers_top_persite_chr6pbs1_25.bed -b annotation_grch38.bed |cut -f1-4,8>HAE_JEJU_pbs_anno25_chr6.bed
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
