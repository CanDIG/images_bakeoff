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

    echo "bwa idx"
    perf-report -o "bwa-${count}" /usr/bin/bwa index hg38.fa.gz

    echo "samtools view"
    perf-report -o "samtools-view-${count}" /usr/bin/samtools view -hf 0x2 "${bam}"  > input.bam

    echo "picard samtofastq"
    perf-report -o "samtofastq-view-${count}" /usr/bin/java -Xmx4g -jar ${PATH_TO_PICARD}/picard.jar SamToFastq VALIDATION_STRINGENCY=LENIENT I=input.bam F=file_1.fastq F2="file_2.fastq"

    echo "bwa mem"
    perf-report -o "bwamem-${count}" /usr/bin/bwa mem hg38.fa.gz file_1.fastq file_2.fastq | samtools view -bS - > output.bam
    
    echo "samtools sort"
    perf-report -o "samtools-sort-${count}" /usr/bin/samtools sort -m 1M output.bam sorted

    rm -- *.fastq output.bam sorted.bam input.bam hg38.fa.*
    count=$((count+1))
done
