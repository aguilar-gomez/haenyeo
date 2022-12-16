#Add names of positions k6
cat header_pos haenyeo.forohanav14.map > positions_selscan
paste positions_selscan scan.haenyeok3.txt > selscan_k3_header.txt

#Filter by likelihood ratio k6
awk '$6 > 2' selscan_k3_header.txt > smallscan_jejuComponent_k3.txt


#Supervised k7 
paste positions_selscan supervisedk7_scan.haenyeok3.txt > supervisedk7_selscan_k3_header.txt
#Filter by likelihood ratio k7
awk '$6 > 2' supervisedk7_selscan_k3_header.txt  > small_supervisedk7_scan.haenyeok3.txt 
