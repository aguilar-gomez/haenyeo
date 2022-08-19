#qpas (population structure)

Kmin=3 # min value of k
Kmax=9 # maximum value of k
mi=450 # max number of iterations
e=0.08 # minimum delta for likelihood

## run ohana structure (qpas)

for k in $(seq $Kmin $Kmax);
do 
echo "running qpas with k $k"
qpas haenyeov113percent.dgm -k $k -qo heanyeo_v11_k${k}_e${e}_mi${mi}_q.matrix -fo heanyeo_v11_k${k}_e${e}_mi${mi}_f.matrix -e ${e} -mi $mi > out.qpas${k}mi${mi}e${e} &
done
