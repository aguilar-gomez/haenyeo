#Convert to CLUES format
chr=$1
pos=$2
grep "\b${pos}\b" ../jeju$chr.haps> genotype.chr$chr.$pos
cut -f6- genotype.chr$chr.$pos -d " "|sed  -e 's/ /\n/g' > chr$chr.$pos.derived
#First step convert data
PATHCLUES=~/programs/CLUES2
python $PATHCLUES/RelateToCLUES.py --RelateSamples chr$chr.$pos.newick --out clues.chr$chr.$pos --DerivedFile chr$
chr.$pos.derived
