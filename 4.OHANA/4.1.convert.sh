OHANA=$PATHTOOHANA

#Recode with plink 
plink --bfile allchr.version14 --recode 12 tab --out haenyeo.forohanav14

#Convert plink format to OHANA
$OHANA/convert ped2dgm haenyeo.forohanav14.ped haenyeov14.dgm 

#Sample sites for OHANA qpas (in this case 3%)
$OHANA/sample-sites.py haenyeov14.dgm  3 ./haenyeov143percent.dgm






