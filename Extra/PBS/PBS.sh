#Max missing genotypes set to 30 (this is version12)
vcftools --recode --vcf selection_scansv11.vcf --max-missing-count 30 --out selection_scansv11_missingCount
Parameters as interpreted:
    --vcf selection_scansv11.vcf
    --max-missing-count 30
    --out selection_scansv11_missingCount
    --recode
#After filtering, kept 5,006,988 out of a possible 6,377,716 Sites

#Subset samples
grep "HAE\|JEJ\|SEO\|CHB\|KOR" threeletterpop.fam > samples4scans
grep "HAE"  samples4scans |cut -f2 >hae.list
grep "JEJ"  samples4scans |cut -f2 >jej.list
grep "SEO\|KOR"  samples4scans |cut -f2 >seo.list
grep "CHB"  samples4scans |cut -f2 >chb.list
grep "HAE\|JEJ"  samples4scans |cut -f2 >jisland.list
cat hae.list jeju.list >haejeju.list

#VCFtools (0.1.17)
#PBS by windows
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

#By site
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
