#Download cramfiles from CHS outgroup
cut -d " " -f2 selscan3popsCHSv2.fam |grep "HG"> CHSindividuals

while read -r ind
do
echo $ind
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CHS/$ind/alignment/$ind.alt_bwamem_GRCh38DH.20150718.CHS.low_coverage.cram{,.crai};
done < CHSindividuals


for file in *.cram; do
  mv "$file" "${file%.alt_bwamem_GRCh38DH.20150718.CHS.low_coverage.cram}.cram"
done


ls /space/s2/diana/korea/haenyeo/GenotypeLikelihoods2023/CHSdata/cramfiles/*cram > CHS.cramlist
ls /space/s2/diana/korea/haenyeo/filtered_fastq/bams38_renamed/*bam > Korean.whole.bamlist

#Remove related individuals:
grep -v "43_S74_L001_chr.bam\|63_S99_L002_chr.bam\|52_S83_L002_chr.bam" Korean.whole.bamlist > unrelated.bamlist

for i in {1..21}; do
  echo "chr$i" >> autosomes
done


cat unrelated.bamlist CHS.cramlist > samples.list4angsd

bamlist=samples.list4angsd
REF=/space/s2/diana/korea/reference/GRCh38/grch38.fasta

/home/diana/bin/angsd -bam $bamlist -out HaeKorChs -minInd 100 -setMinDepthInd 1 -minMapQ 25 \
-minQ 25 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -GL 1 -doMaf 1 -doMajorMinor 4 -doGlf 2 \
-SNP_pval 1e-6 -nThreads 16 -ref $REF -setMaxDepthInd 50 -doCounts 1 -skipTriallelic 1 -dosnpstat 1 \
-doHWE 1 -sb_pval 1e-4 -hetbias_pval 1e-6 -rf autosomes \
-edge_pval 1e-4 -mapQ_pval 1e-4
