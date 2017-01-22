#!/usr/bin/env bash

readonly DATA="${HOME}/data"
readonly PATH_TO_PICARD="/usr/local/share/java"

let count=1
for bam in ${DATA}/*.bam
do
    echo "$count: $bam"
    outdir="${HOME}/out/${count}"
    mkdir -p "${outdir}"
    cd "${outdir}" || exit

    cp "${HOME}/data/hg38.fa.gz" .
    cp "${HOME}/data/hg38.fa.gz.fai" .

    if [[ "$bam" =~ chrom20 ]]
    then
        viewscript=${DATA}/view_exons_20.sh
    else
        viewscript=${DATA}/view_exons.sh
    fi
    cp ${viewscript} ./samtools_view.sh

    echo "bwa idx"
    time bwa index hg38.fa.gz 

    echo "samtools view"
    time samtools view -hf 0x2 ${bam}  > input.bam

    echo "picard samtofastq"
    time java -Xmx4g -jar ${PATH_TO_PICARD}/picard.jar SamToFastq VALIDATION_STRINGENCY=LENIENT I=input.bam F=file_1.fastq F2="file_2.fastq"

    echo "bwa mem"
    time bwa mem hg38.fa.gz file_1.fastq file_2.fastq | samtools view -bS - > output.bam
    
    echo "samtools sort"
    time samtools sort -m 1M output.bam sorted

    echo "samtools index"
    time samtools index sorted.bam

    echo "samtools view random"
    time ./samtools_view.sh 

    echo "rm"
    time rm *.fastq output.bam sorted.bam* input.bam hg38.fa.* >&2 rm-${count}.txt
    count=$((count+1))
done
