OHANA=$PATHTOOHANA

#Recode with plink 
plink --bfile allchr.version13 --recode 12 tab --out haenyeo.forohanav13 

#Convert plink format to OHANA
$OHANA/convert ped2dgm haenyeo.forohanav11.ped haenyeov11.dgm 
$OHANA/convert ped2dgm haenyeo.forohanav13.ped haenyeov13.dgm 

#Sample sites for OHANA qpas (in this case 3%)
$OHANA/sample-sites.py haenyeov13.dgm  3 ./haenyeov133percent.dgm






