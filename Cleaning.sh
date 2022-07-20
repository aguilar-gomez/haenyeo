#Generate fastqc
for fastqfile in *fastq.gz
do 
fastqc -o /space/s2/diana/haenyeo/data/fastqc -f fastq -t 5 $fastqfile 
done

#Clean reads PRINSEQ-lite 0.20.4
for sample in *R1*fastq
do  name=${sample%%_R*}
    echo $name
            prinseq-lite.pl -fastq ${name}"_R1"* -fastq2 ${name}"_R2"* -lc_method dust -lc_threshold 7 -custom_params "G 50" -out_good $name.filter
            fastqc -o /space/s2/diana/haenyeo/data/fastqc -f fastq -t 10 $name.filter_1.fastq
            fastqc -o /space/s2/diana/haenyeo/data/fastqc -f fastq -t 10 $name.filter_2.fastq
    gzip ${name}"_R1"* ${name}"_R2"*
            rm *bad*
            rm *single*
done
