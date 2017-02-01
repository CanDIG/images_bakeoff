#!/usr/bin/env bash

readonly DATA="${HOME}/data"
readonly OUT="${HOME}/singularity_out"

function log {
    local file=$1
    local msg=$2
    echo "$msg" 
    echo "$msg" >> "$file"
}

readonly BWA=${HOME}/singularities/ubuntu14-bwa.img
readonly SAMTOOLS=${HOME}/singularities/ubuntu14-samtools-1.3.1.img
readonly PICARD=${HOME}/singularities/ubuntu16-picard.img

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
    /usr/bin/time -p -a -o "$logfile" singularity run ${BWA} index hg38.fa.gz 

    log "$logfile" "samtools view"
    /usr/bin/time -p -a -o "$logfile" singularity run ${SAMTOOLS} view -hf 0x2 /data/$( basename ${bam} )  > input.bam

    log "$logfile" "picard samtofastq"
    /usr/bin/time -p -a -o "$logfile" singularity run ${PICARD} VALIDATION_STRINGENCY=LENIENT I=input.bam F=file_1.fastq F2="file_2.fastq"

    log "$logfile" "bwa mem"
    /usr/bin/time -p -a -o "$logfile" singularity run ${BWA} mem hg38.fa.gz file_1.fastq file_2.fastq | singularity run ${SAMTOOLS} view -bS - > output.bam
    
    log "$logfile" "samtools sort"
    /usr/bin/time -p -a -o "$logfile" singularity run ${SAMTOOLS} sort -m 50M output.bam -o sorted.bam

    log "$logfile" "samtools index"
    /usr/bin/time -p -a -o "$logfile" singularity run ${SAMTOOLS} index sorted.bam

    log "$logfile" "samtools view random"
    /usr/bin/time -p -a -o "$logfile" singularity exec ${SAMTOOLS} /usr/local/bin/samtools_random.sh viewlines.txt | wc -l

    log "$logfile" "rm"
    /usr/bin/time -p -a -o "$logfile" singularity exec ${BWA} /bin/rm -ubuntu -f file_1.fastq file_2.fastq output.bam sorted.bam sorted.bam.bai input.bam hg38.fa.gz hg38.fa.gz.amb hg38.fa.gz.ann hg38.fa.gz.bwt hg38.fa.gz.fai hg38.fa.gz.pac hg38.fa.gz.sa viewlines.txt

    count=$((count+1))
done
