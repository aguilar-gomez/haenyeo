460 individuals try for k=6 and k=7
sed  's/HAE/0/' PopIndv>temp
sed -i 's/JEJ/1/' temp
sed -i 's/SEO/2/' temp
sed -i 's/GBR/3/' temp
sed -i 's/KHV/4/' temp
sed -i 's/KOR/5/' temp
sed -i 's/CHS/6/' temp
sed -i 's/CHB/7/' temp
sed -i 's/JPT/8/' temp
sed -i 's/YRI/9/' temp
tr -d "\n\r" < temp | sed 's/./& /g' | fold -w80 >popasign.mat


k=7 
mi=450
e=0.08
qpas haenyeov143percent.dgm -k $k -qo supervised_v14_k${k}_e${e}_mi${mi}_q.matrix -fo supervised_v14_k${k}_e${e}_mi${mi}_f.matrix -e ${e} -mi $mi -fg fgk${k} > supervised.out.qpas${k}mi${mi}e${e} &
