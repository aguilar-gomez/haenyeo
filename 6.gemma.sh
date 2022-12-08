#Extract Samples collected by Melissa
 plink --file haenyeo.forohanav14 --out assotest_haenyeov14 --keep samples2keep --make-bed

#Modify the fam file to include phenotypes
cp /space/s2/diana/korea/haenyeo/June2022/plink/gemma/assotest_haenyeov12.fam assotest_haenyeov14.fam

#Calculate relatedness matrix with all positions
gemma -bfile assotest_haenyeov14 -gk 2 -o relatedness_matrix.haenyeo

#Extract positions
 plink --bfile assotest_haenyeov14 --out assotest_haenyeov14_pos2test --extract SNP2test --make-bed

#Get covariate file
cp /space/s2/diana/korea/haenyeo/June2022/plink/gemma/diving_Cov .

#Run association
gemma -bfile assotest_haenyeov14_pos2test -lmm 2 -o lrt_top10SNPs_diveCov -c diving_Cov -k output/relatedness_matrix.haenyeo.sXX.txt

gemma -bfile assotest_haenyeov14_pos2test -lmm 2 -o lrt_top10SNPs_diveCov_forceSNPs -c diving_Cov -k output/relatedness_matrix.haenyeo.sXX.txt -notsnp -miss 1
