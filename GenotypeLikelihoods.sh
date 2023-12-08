#Download cramfiles from CHS outgroup
cut -d " " -f2 selscan3popsCHSv2.fam |grep "HG"> CHSindividuals

while read -r ind
do
echo $ind
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CHS/$ind/alignment/$ind.alt_bwamem_GRCh38DH.20150718.CHS.low_coverage.cram{,.crai};
done < CHSindividuals



