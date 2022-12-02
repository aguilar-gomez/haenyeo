Kmin=4 # min value of k
Kmax=8 # maximum value of k
mi=450 # max number of iterations
e=0.08 # minimum delta for likelihood

for k in $(seq $Kmin $Kmax);
do
echo "running nemeco with k $k"
$OHANA/nemeco haenyeov143percent.dgm heanyeo_v14_k${k}_e${e}_mi${mi}_f.matrix -co c.matrix_k${k} -mi 5 > out1.nemk${k} 
$OHANA/convert cov2nwk c.matrix_k${k} haenyeov14_k${k}.nwk
tail -n +2 heanyeo_v14_k${k}_e${e}_mi${mi}_q.matrix  > haenyeov14_k${k}.Q
done


#supervised
 tail -n +2 supervised_v14_k7_e0.08_mi450_q.matrix  > haenyeov14_sup_k${k}.Q
