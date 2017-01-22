 #!/bin/bash

function usage 
{
    >&2 echo "$0 /path/to/data/dir"
    >&2 echo "   downloads needed data for benchmark."
    exit 1
}

function download_exons {
    local datadir=$1

    if [[ ! -f ${datadir}/exons.bed ]]
    then
         curl  -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz" \
             | gunzip -c \
             | awk '{n=int($8); split($9,S,/,/);split($10,E,/,/); for(i=1;i<=n;++i) {printf("%s\t%s\t%s\n",$2,S[i],E[i]);} }' \
             | sed -e 's/chr//' \
             > "${datadir}"/exons.bed
    fi

    awk '$1==20' ${datadir}/exons.bed \
        | uniq \
        | sort -R \
        > "${datadir}"/exons_20-random.bed

    sort -R ${datadir}/exons.bed \
        | awk 'BEGIN{ for(i=1;i<=22;i++) {ischr[i ""]="1"}; ischr["X"]="1"; ischr["Y"]="1";} ischr[$1 ""]=="1"' \
        | uniq \
        | sort -R \
        > "${datadir}"/exons-random.bed
}

function output_view {
    local input=$1
    local output=$2
    local perline=${3:-1000}

    echo "#!/usr/bin/env bash" > ${output}
    nentries=$( wc -l "$input" | cut -f 1 -d ' ')
    nlines=$((nentries/perline))

    for i in $( seq 1 "${nlines}" )
    do
        last=$((i*perline))
        records=$( cat $input | head -n $last | tail -n $nentries | awk '{printf "%s:%s-%s ", $1, $2, $3}' )
        echo "/usr/bin/samtools view -F 4 sorted.bam ${records} | wc -l" >> "${output}"
    done
}

function main {
    local datadir=$1
    local nperline=${2:-1000}
    download_exons "$datadir"
    output_view "$datadir"/exons_20-random.bed "$datadir"/view_exons_20.sh $nperline
    output_view "$datadir"/exons-random.bed "$datadir"/view_exons.sh $nperline
}

if [[ -z "$1" || ! -d "$1" ]]
then
    usage
fi
main $1 1000
