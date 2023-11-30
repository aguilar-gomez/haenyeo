
#k6 
#Generate matrices with value added
for k in $(seq 0 5);
do
cmatrixK.py c.matrix_k6 $k 10
done
 
#Run qpas to calculate allele frequency matrix for all positions using the inferred components with k=6
$OHANA/qpas haenyeov14.dgm -k 6 -qi heanyeo_v14_k6_e0.08_mi450_q.matrix -fo haenyeov14_k6.full.f.matrix -e 0.08 -mi 450 >out.qpas.full &

#Run selscan
for k in $(seq 0 5);
do 
echo "running selscan with selection in component k $k"
$OHANA/selscan haenyeov14.dgm haenyeov14_k6.full.f.matrix c.matrix_k6 -cs c.matrix_k6selK${k}_h10 > scan.haenyeok$k.txt &
done
