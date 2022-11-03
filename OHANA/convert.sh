OHANA=$PATHTOOHANA

#Recode with plink 
plink --bfile pcapops.allchr.maf2.v11 --recode 12 tab --out haenyeo.forohanav11 

#Convert plink format to OHANA
$OHANA/convert ped2dgm haenyeo.forohanav11.ped haenyeov11.dgm 

#Sample sites for OHANA qpas (in this case 3%)
$OHANA/sample-sites.py haenyeov11.dgm  3 ./haenyeov113percent.dgm






