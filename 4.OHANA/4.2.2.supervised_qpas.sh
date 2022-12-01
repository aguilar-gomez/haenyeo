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
