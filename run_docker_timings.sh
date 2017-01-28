#!/usr/bin/env bash

readonly DATA="${HOME}/data"
readonly OUT="${HOME}/docker_out"

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
        viewscript=${DATA}/records_20.txt
    else
        viewscript=${DATA}/records.txt
    fi
    cp ${viewscript} ./viewlines.txt

    log "$logfile" "bwa idx"
    /usr/bin/time -p -a -o "$logfile" docker run -v ${PWD}:/pwd bwa-ubuntu index /pwd/hg38.fa.gz 

    log "$logfile" "samtools view"
    /usr/bin/time -p -a -o "$logfile" docker run -v ${PWD}:/pwd -v ${DATA}:/data samtools-alpine view -hf 0x2 /data/$( basename ${bam} )  > input.bam

    log "$logfile" "picard samtofastq"
    /usr/bin/time -p -a -o "$logfile" docker run -v ${PWD}:/pwd picard SamToFastq VALIDATION_STRINGENCY=LENIENT I=/pwd/input.bam F=/pwd/file_1.fastq F2="/pwd/file_2.fastq"

    log "$logfile" "bwa mem"
    /usr/bin/time -p -a -o "$logfile" docker run -v ${PWD}:/pwd bwa-ubuntu mem /pwd/hg38.fa.gz /pwd/file_1.fastq /pwd/file_2.fastq > output.sam && docker run -v ${PWD}:/pwd samtools-alpine view -bS /pwd/output.sam > output.bam && rm output.sam
    
    log "$logfile" "samtools sort"
    /usr/bin/time -p -a -o "$logfile" docker run -v ${PWD}:/pwd samtools-alpine sort -m 50M /pwd/output.bam -o /pwd/sorted.bam

    log "$logfile" "samtools index"
    /usr/bin/time -p -a -o "$logfile" docker run -v ${PWD}:/pwd samtools-alpine index /pwd/sorted.bam

    log "$logfile" "samtools view random"
    /usr/bin/time -p -a -o "$logfile" docker run -v ${PWD}:/pwd --entrypoint /usr/wrapper/samtools_random samtools-alpine /pwd/sorted.bam /pwd/viewlines.txt | wc -l

    log "$logfile" "rm"
    /usr/bin/time -p -a -o "$logfile" docker run -v ${PWD}:/pwd --entrypoint /bin/rm bwa-ubuntu -f /pwd/{file_1.fastq,file_2.fastq,output.bam,sorted.bam,sorted.bam.bai,input.bam,hg38.fa.gz,hg38.fa.gz.amb,hg38.fa.gz.ann,hg38.fa.gz.bwt,hg38.fa.gz.fai,hg38.fa.gz.pac,hg38.fa.gz.sa,viewlines.txt}

    count=$((count+1))
done
