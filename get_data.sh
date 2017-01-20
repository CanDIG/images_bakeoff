#!/usr/bin/env bash

function usage 
{
    >&2 echo "$0 /path/to/data/dir"
    >&2 echo "   downloads needed data for benchmark."
    exit 1
}

readonly BAMFILES=(ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00119/alignment/HG00119.chrom20.ILLUMINA.bwa.GBR.low_coverage.20101123.bam \
                   ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/NA19700/alignment/NA19700.chrom20.ILLUMINA.bwa.ASW.low_coverage.20101123.bam \
                   ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00513/alignment/HG00513.chrom20.ILLUMINA.bwa.CHS.low_coverage.20101123.bam \
                   ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00119/alignment/HG00119.unmapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam \
                   ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/NA19700/alignment/NA19700.unmapped.ILLUMINA.bwa.ASW.low_coverage.20101123.bam \
                   ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00513/alignment/HG00513.unmapped.ILLUMINA.bwa.CHS.low_coverage.20101123.bam \
                   ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00119/alignment/HG00119.mapped.ILLUMINA.bwa.GBR.low_coverage.20101123.bam \
                   ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/NA19700/alignment/NA19700.mapped.ILLUMINA.bwa.ASW.low_coverage.20101123.bam \
                   ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00513/alignment/HG00513.unmapped.ILLUMINA.bwa.CHS.low_coverage.20101123.bam \
                   ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/data/CHS/HG00513/high_cov_alignment/HG00513.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.cram)

readonly REFERENCE=http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

function executable_exists 
{
    local exe=$1
    hash "$exe" 2>/dev/null
}

function gnuunzipcat {
    if executable_exists "gzcat"
    then
        gzcat "$@"
    else
        zcat "$@"
    fi
}

function download_file 
{
    local url=$1
    local default_output=$( basename "$url" )
    local output=${2:-${default_output}}
    local NTHREADS=4

    # download using axel, wget, or curl, in order of preference
    echo "Downloading $url to $output"
    if executable_exists "axel"
    then
        axel -q -n "$NTHREADS" "$url" -o "$output"
    elif executable_exists "wget"
    then
        wget -nv "$url" -O "$output"
    elif executable_exists "curl"
    then
        curl -S "$url" -o "$output"
    else
        >&2 echo "$0: Neither axel, wget, nor curl found in PATH ${PATH}"
        exit 1
    fi
}

function download_ref 
{
    local url=$1
    local data_path=$2
    local output=${3:-$( basename "$url" )}

    output_path="${data_path}/${output}"
    tmp_path="${data_path}/tmp_${output}"

    if ! executable_exists "bgzip"  || ! executable_exists "samtools"
    then
        >&2 echo "$0: need bgzip and samtools to process reference"
        exit 1
    fi

    download_file "$url" "$tmp_path"
    gnuunzipcat "$tmp_path" \
        | bgzip \
        > "${output_path}" \
        && rm "$tmp_path" 

    echo "Indexing reference..."
    samtools faidx "${output_path}"
}

function main 
{
    local data_path="$1"

    if [[ ! -d "$data_path" ]] 
    then
        >&2 echo "Invalid path: $data_path"
        exit
    fi
    cd "$data_path" || exit 1
    download_ref "$REFERENCE" "." hg38.fa.gz
    for bam in "${BAMFILES[@]}"
    do
        download_file "$bam"
    done
}

if [[ -z "$1" ]] || [[ ! -d "$1" ]]
then
	usage
fi
main "$1"
