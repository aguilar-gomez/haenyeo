
relatemac/bin/RelateFileFormats --mode ConvertFromVcf --haps Step1.haps --sample Step1.sample -i chr1.subsetKorean
##Here, snps with any missing data are automatically filtered out, this filters out approximately 1/2 of all the snps, possibly including SNPs covered by the mask.

##Ancestral sequence comes from: https://ftp.ensembl.org/pub/release-109/fasta/ancestral_alleles/
##Masks come from: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20160622_genome_mask_GRCh38/StrictMask/

relatemac/scripts/PrepareInputFiles/PrepareInputFiles.sh  --haps Step1.haps  --sample Step1.sample  --ancestor homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_1.fa --mask  20160622.chr1.mask.fasta -o Step2 --poplabels allSamplesOrder.popfile
#In this step, we go from 1297478 SNPs to 859977 SNPs


#Genetic map comes from: https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/ then, we split on chromosomes


time relatemac/bin/Relate  --mode All   -m 1.25e-8     -N 30000  \
      --haps Step2.haps.gz \
      --sample Step2.sample.gz \
      --map maps/map1.txt  \
      --seed 5 \
      -o Chrom1 \
      --annot Step2.annot \
      --dist Step2.dist.gz 

relatemac/scripts/EstimatePopulationSize/EstimatePopulationSize.sh  -i Chrom1 -m 1.25e-8 --years_per_gen 27 --threads 4 --seed 43 --output Chr1Inferred_27years --poplabels allSamplesOrder.popfile  --pop_of_interest KOR,JEJ

rm *_avg.rate
rm *.pairwise.bin

relatemac/bin/RelateExtract --mode SubTreesForSubpopulation  --anc Chr1Inferred_27years.anc.gz  --mut Chr1Inferred_27years.mut.gz --poplabels allSamplesOrder.popfile   --pop_of_interest JEJ --output Chr1_JEJ
relatemac/bin/RelateExtract --mode SubTreesForSubpopulation  --anc Chr1Inferred_27years.anc.gz  --mut Chr1Inferred_27years.mut.gz --poplabels allSamplesOrder.popfile   --pop_of_interest KOR --output Chr1_KOR

##Here, we have to curate the coal files by hand into KOR.coal and JEJ.coal

baseposition=161704797 #the relevant snp
numsampss=200 #the number of samples we take

relatemac/scripts/SampleBranchLengths/SampleBranchLengths.sh  --mu 1.25e-8 \
                 -i Chr1_JEJ \
                 -o JEJ_Samples \
                 --coal JEJ.coal \
                 --format n \
                 --num_samples  $numsampss \
                 --first_bp  $baseposition \
                 --last_bp  $baseposition \
                 --dist Chr1Inferred_27years.dist \
                 --seed 12

relatemac/scripts/SampleBranchLengths/SampleBranchLengths.sh  --mu 1.25e-8 \
                 -i Chr1_KOR \
                 -o KOR_Samples \
                 --coal KOR.coal \
                 --format n \
                 --num_samples  $numsampss \
                 --first_bp  $baseposition \
                 --last_bp  $baseposition \
                 --dist Chr1Inferred_27years.dist \
                 --seed 12


grep ${baseposition} Step1.haps > derivedfile.txt
sed -i -e 's/ /\n/g' derivedfile.txt
sed -i -e 1,5d derivedfile.txt
rm derivedfile.txt-e
head -n 102  derivedfile.txt > JEJDerivedFile.txt # because first 51 diploid samples are jeju
tail -n 216  derivedfile.txt > KORDerivedFile.txt

python3.9 CLUES2/RelateToCLUES.py  --RelateSamples KOR_Samples.newick --DerivedFile KORDerivedFile.txt --out KOR_Times.txt # this flips 2 out of 216 leaves
python3.9 CLUES2/RelateToCLUES.py  --RelateSamples JEJ_Samples.newick --DerivedFile JEJDerivedFile.txt --out JEJ_Times.txt #this flips 7 out of 102 leaves


KORderivedfreq=0.0787037
JEJderivedfreq=0.3333333333
tcutofff=450

python3.10 ~/desktop/CLUES2/inference.py --coal KOR.coal --popFreq  ${KORderivedfreq} --times KOR_Times.txt_times.txt  --out KOR_results  --tCutoff ${tcutofff}  --noAlleleTraj
python3.10 ~/desktop/CLUES2/inference.py --coal JEJ.coal --popFreq  ${JEJderivedfreq} --times JEJ_Times.txt_times.txt  --out JEJ_results  --tCutoff ${tcutofff}  --noAlleleTraj

python3.10 ~/desktop/CLUES2/inference.py --coal JEJ.coal --popFreq  ${JEJderivedfreq} --times JEJ_Times.txt_times.txt --out JEJ_results_1215  --tCutoff ${tcutofff} --timeBins 45 #1215 years
python3.10 ~/desktop/CLUES2/inference.py --coal JEJ.coal --popFreq  ${JEJderivedfreq} --times JEJ_Times.txt_times.txt --out JEJ_results_6993  --tCutoff ${tcutofff}  --noAlleleTraj --timeBins 259 #6993 years ago.

python3.10 ~/desktop/CLUES2/plot_traj.py --freqs  JEJ_results_1215_freqs.txt --post JEJ_results_1215_post.txt --figure JEJ_results_1215
