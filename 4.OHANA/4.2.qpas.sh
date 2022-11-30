#qpas (population structure)

Kmin=4 # min value of k
Kmax=8 # maximum value of k
mi=450 # max number of iterations
e=0.08 # minimum delta for likelihood

## run ohana structure (qpas)

for k in $(seq $Kmin $Kmax);
do 
echo "running qpas with k $k"
$OHANA/qpas haenyeov143percent.dgm -k $k -qo heanyeo_v14_k${k}_e${e}_mi${mi}_q.matrix -fo heanyeo_v14_k${k}_e${e}_mi${mi}_f.matrix -e ${e} -mi $mi > out.qpas${k}mi${mi}e${e} &
done
