#Generate matrices with value added
for k in $(seq 0 7);
do
cmatrixK.py c.matrix_k8 $k 10
done
 
#Run qpas to calculate allele frequency matrix for all positions using the inferred components with k=7
$OHANA/qpas haenyeov13.dgm -k 7 -qi heanyeo_v13_k8_e0.08_mi450_q.matrix -fo haenyeov13_k7.full.f.matrix -e 0.08 -mi 450 >out.qpas.full &

#Run selscan
for k in $(seq 0 6);
do 
echo "running selscan with selection in component k $k"
$OHANA/selscan haenyeov13.dgm haenyeov13_k7.full.f.matrix c.matrix_k7 -cs c.matrix_k7selK${k}_h10 > scan.haenyeok$k.txt &
done
