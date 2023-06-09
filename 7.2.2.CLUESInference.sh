chr=$1
pos=$2
cutoff=$3
average=$(awk '{ sum += $1 } END { print sum / NR }' chr$chr.$pos.derived)

echo Frequency of derived allele $average

PATHCLUES=~/programs/CLUES2
python $PATHCLUES/inference.py \
--times clues.chr$chr.${pos}_times.txt \
--popFreq $average \
--coal ../jeju.$chr.popsize.coal \
--tCutoff $cutoff \
--df 400 \
--out chr$chr.$pos.t$cutoff



