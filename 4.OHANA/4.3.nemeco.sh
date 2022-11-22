Kmin=3 # min value of k
Kmax=9 # maximum value of k
mi=450 # max number of iterations
e=0.08 # minimum delta for likelihood

for k in $(seq $Kmin $Kmax);
do
echo "running nemeco with k $k"
$OHANA/nemeco haenyeov133percent.dgm heanyeo_v13_k${k}_e${e}_mi${mi}_f.matrix -co c.matrix_k${k} -mi 5 > out1.nemk${k} 
$OHANA/convert cov2nwk c.matrix_k${k} haenyeov13_k${k}.nwk
tail -n +2 heanyeo_v13_k${k}_e${e}_mi${mi}_q.matrix  > haenyeov13_k${k}.Q
done
