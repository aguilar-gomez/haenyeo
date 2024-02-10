FILE=haeProject
for i in {0..5}
do
 treemix -i $FILE.treemix.frq.gz -m $i -o $FILE.$i -bootstrap -k 500 -root GBR > treemix_${i}_log &
done
