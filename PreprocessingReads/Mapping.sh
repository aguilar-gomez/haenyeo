for sample in *.filter_1*
do
        name=${sample%.filter*}
        name_o=${sample%_*}
        echo $name_o"_1"*
        echo $name_o"_2"*
        bwa mem -t 8 /media/external3/diana/reference/GRCh38/GRCh38_genomic.fna $name_o"_1"* $name_o"_2"* >$name.sam
        samtools view -Sb -F 1804 $name.sam|samtools sort -@ 20 - >$name.s38.bam
        samtools index $name.s38.bam
        rm $name.sam
done
 
 
Coverage/depth
for bam in *bam
do
  samtools depth -b /media/external3/diana/reference/chr_autosomes.bed -o $bam.depth.autosomes $bam
Done
 
 
for bam in *bam
do
echo $bam
  samtools coverage -o $bam.cov $bam
done
