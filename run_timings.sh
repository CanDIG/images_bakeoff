#!/usr/bin/env bash

readonly DATA="${HOME}/data"
readonly OUT="${HOME}/out"
readonly PATH_TO_PICARD="/usr/local/share/java"

function log {
    local file=$1
    local msg=$2
    echo "$msg" 
    echo "$msg" >> "$file"
}

let count=1
for bam in ${DATA}/*.bam
do
    echo "$count: $bam"
    outdir="${OUT}/${count}"
    mkdir -p "${outdir}"
    cd "${outdir}" || exit

    logfile="./timings.log"

    cp "${DATA}/hg38.fa.gz" .
    cp "${DATA}/hg38.fa.gz.fai" .
    if [[ "$bam" =~ chrom20 ]]
    then
        viewscript=${DATA}/view_exons_20.sh
    else
        viewscript=${DATA}/view_exons.sh
    fi
    cp "${viewscript}" ./samtools_view.sh
    chmod 755 ./samtools_view.sh

    log "$logfile" "bwa idx"
    /usr/bin/time -p -a -o "$logfile" bwa index hg38.fa.gz 

    log "$logfile" "samtools view"
    /usr/bin/time -p -a -o "$logfile" samtools view -hf 0x2 "${bam}"  > input.bam

    log "$logfile" "picard samtofastq"
    /usr/bin/time -p -a -o "$logfile" java -Xmx4g -jar ${PATH_TO_PICARD}/picard.jar SamToFastq VALIDATION_STRINGENCY=LENIENT I=input.bam F=file_1.fastq F2="file_2.fastq"

    log "$logfile" "bwa mem"
    /usr/bin/time -p -a -o "$logfile" bwa mem hg38.fa.gz file_1.fastq file_2.fastq | samtools view -bS - > output.bam
    
    log "$logfile" "samtools sort"
    /usr/bin/time -p -a -o "$logfile" samtools sort -m 50M output.bam sorted

    log "$logfile" "samtools index"
    /usr/bin/time -p -a -o "$logfile" samtools index sorted.bam

    log "$logfile" "samtools view random"
    /usr/bin/time -p -a -o "$logfile" ./samtools_view.sh 

    log "$logfile" "rm"
    /usr/bin/time -p -a -o "$logfile" rm -f -- *.fastq output.bam sorted.bam* input.bam hg38.fa.* 
    count=$((count+1))
done
