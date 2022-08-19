Kmin=3 # min value of k
Kmax=9 # maximum value of k

for k in $(seq $Kmin $Kmax);
do
echo "running nemeco with k $k"
$OHANA/nemeco haenyeov113percent.dgm heanyeo_v11_k${k}_e${e}_mi${mi}_f.matrix -co c.matrix_k${k} -mi 5 > out1.nemk${k} 
$OHANA/convert cov2nwk c.matrix_k${k} haenyeov11_k${k}.nwk
tail -n +2 heanyeo_v11_k${k}_e${e}_mi${mi}_q.matrix  > haenyeov11_k${k}.Q
done
