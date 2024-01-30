#Extract samples from files with all data, no SNP filtering
for chr in chr{1..22}
do
bcftools view -S allKorean $chr.allK2023.vcf.gz >$chr.subsetKorean.vcf
done

nohup bash extractVCF.sh > out.extractVCF &

#Convert from vcf 

#Convert to Relate Format

PATH_TO_RELATE=~/programs/relate_ARGS/
for chr in chr{1..22}
do
$PATH_TO_RELATE/bin/RelateFileFormats \
	--mode ConvertFromVcf \ 	
	--haps $chr.haps  	\
	--sample $chr.sample \	
	-i  $chr.subsetKorean 
done

nohup bash convertVCF.sh > out.covertVCF &


#Obtain trees
PATH_TO_RELATE=~/programs/relate_ARGS/
for chr in chr{1..22}
do
$PATH_TO_RELATE/bin/Relate \
      --mode All \
      -m 1.25e-8 \
      -N 10000 \
      --haps $chr.haps \
      --sample $chr.sample \
      --map maps/map$chr.txt \
      --seed 1 \
      -o relate_${chr} &
done

nohup bash ObtainTrees.sh > out.ObtainTrees &


#Estimate popSize

PATH_TO_RELATE=~/programs/relate_ARGS/
$PATH_TO_RELATE/scripts/EstimatePopulationSize/EstimatePopulationSize.sh \
              -i relate \
              --threads 10 \
              -m 1.25e-8 \
              --years_per_gen 27 \
              --seed 1 \
              --poplabels allSamplesOrder.popfile \
              --first_chr 1 \
              --last_chr 22 \
              -o gentime27AllChr

nohup sh EstimatePopSize.sh > out.popsize &

#Extract Subtree
PATH_TO_RELATE=~/programs/relate_ARGS/
$PATH_TO_RELATE/bin/RelateExtract\
                 --mode SubTreesForSubpopulation\
                 --anc relate_chr1.anc \
                 --mut relate_chr1.mut \
                 --poplabels allSamplesOrder.popfile \
                 --pop_of_interest JEJ \
                 -o JEJ_chr1

#Sample Branch Length

chr=$1
pos=$2

PATH_TO_RELATE=~/programs/relate_ARGS/
$PATH_TO_RELATE/scripts/SampleBranchLengths/SampleBranchLengths.sh \
                -i JEJ_chr1 \
                -o chr$chr.$pos.JEJ  \
                -m 1.25e-8 \
                --coal gentime27AllChr.coal \
                --format n \
                --num_samples 5 \
                --first_bp $pos \
                --last_bp $pos

nohup bash SampleBL.sh 1 161704797 > out.snp1 &


pos=161704797 
chr=1
 
grep $pos chr$chr.haps> genotype.chr$chr.$pos
cut -f6- genotype.chr$chr.$pos -d " "|sed  -e 's/ /\n/g' > chr$chr.$pos.derived

#The top 51 samples in the file are JEJ
head -102 chr$chr.$pos.derived >chr$chr.$pos.derived.JEJ

#Convert to CLUES format

PATHCLUES=~/programs/CLUES2
python $PATHCLUES/RelateToCLUES.py --RelateSamples chr$chr.$pos.JEJ.newick --out clues.chr$chr.$pos --DerivedFile chr$chr.$pos.derived.JEJ

Second step inference
chr=$1
pos=$2
cutoff=$3
step=$4
average=$(awk '{ sum += $1 } END { print sum / NR }' chr$chr.$pos.derived)
echo Frequency of derived allele $average
PATHCLUES=~/programs/CLUES2
python $PATHCLUES/inference.py \
	--times clues.chr$chr.${pos}_times.txt \
	--popFreq $average \
	--coal gentime27AllChr.coal \
	--tCutoff $cutoff \
	--df 400 \
	--out chr$chr.$pos.t$cutoff.bin$step \
	--timeBins $step


nohup sh CLUESInference.sh $chr $pos 1000 100 > out.clues &
nohup sh CLUESInference.sh $chr $pos 1000 200 > out.clues &

PATHCLUES=~/programs/CLUES2

python3  $PATHCLUES/plot_traj.py --freqs chr1.161704797.t1000.bin100_freqs.txt --post chr1.161704797.t1000.bin100_post.txt --figure t1000bin100 --generation_time 27 &

python3  $PATHCLUES/plot_traj.py --freqs chr1.161704797.t1000.bin200_freqs.txt --post chr1.161704797.t1000.bin200_post.txt --figure t1000bin200 --generation_time 27 &

