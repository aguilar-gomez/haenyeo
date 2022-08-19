OHANA=$PATHTOOHANA

#Convert plink format to OHANA
$OHANA/convert ped2dgm haenyeo.forohanav11.ped haenyeov11.dgm 

#Sample sites for OHANA qpas (in this case 3%)
$OHANA/sample-sites.py haenyeov11.dgm  3 ./haenyeov113percent.dgm






