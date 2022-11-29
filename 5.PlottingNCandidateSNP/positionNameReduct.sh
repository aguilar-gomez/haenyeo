#Add names of positions
cat header_pos haenyeo.forohanav13.map > positions_selscan
paste positions_selscan scan.haenyeok4.txt > selscan_k4_header.txt

#Filter by likelihood ratio
awk '$6 > 2' selscan_k4_header.txt > smallscan_jejuComponent_k4.txt
