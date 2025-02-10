#Extract Samples collected by Melissa
 plink --file haenyeo.forohanav14 --out assotest_haenyeov14 --keep samples2keep --make-bed

#Modify the fam file to include phenotypes
cp /space/s2/diana/korea/haenyeo/June2022/plink/gemma/assotest_haenyeov12.fam assotest_haenyeov14.fam

#Calculate relatedness matrix with all positions
gemma -bfile assotest_haenyeov14 -gk 2 -o relatedness_matrix.haenyeo

#Extract positions
plink --bfile assotest_haenyeov14 --out assotest_haenyeov14_pos2test --extract SNP2test --make-bed
plink --bfile assotest_haenyeov14 --out assotest_haenyeov14_allanno --extract posallanno --make-bed

#Run association
plink --bfile asso_haenyeo_pos2test --linear --covar covs_for_plink --out dias_untrans_top10 --allow-no-sex --hide-covar
