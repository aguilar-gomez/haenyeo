#Convert to format to use in Treemix and Admixture Bayes
awk '{ print $1, $2, $1 }' OFS="\t"  allchrfilteredAllMaf1.fam > pops.clust

#Convert SEO,KOA and KMJ to KOR
sed -i 's/SEO/KOR/g' pops.clust
sed -i 's/KOA/KOR/g' pops.clust
sed -i 's/KMJ/KOR/g' pops.clust

#Convert JEU and HAE to JEJ
sed -i 's/JEU/JEJ/g' pops.clust
sed -i 's/HAE/JEJ/g' pops.clust

#Remove YRI samples, because of Neanderthal introgression mess
grep "YRI" allchrfilteredAllMaf1.fam > YRI
plink --bfile allchrfilteredAllMaf1 --geno 0 --out haeProject --double-id --remove-fam YRI --make-bed

#Generate frq file
plink --bfile haeProject --geno 0 --freq --missing --within  pops.clust --out haeProject --double-id --remove-fam YRI 
gzip haeProject.frq.strat
plink2treemix.py haeProject.frq.strat.gz haeProject.treemix.frq.gz
