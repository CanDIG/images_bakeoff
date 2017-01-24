#!/bin/bash
##
## wrapper script for abusing samtools.  
## 

readonly BAMFILE=$1
readonly REGIONSFILE=$2
readonly SAMTOOLS=/usr/local/bin/samtools

if [[ -z "$BAMFILE" ]] || [[ ! -f "$BAMFILE" ]] || [[ -z "$REGIONSFILE" ]] || [[ ! -f "$REGIONSFILE" ]] 
then
    >&2 echo "Usage: $0 reads.bam regions.txt"
    >&2 echo "  repeatedly runs samtools view -F 4 reads.bam regions..."
    >&2 echo "  where each regions is a line from the regions file."
    exit 1
fi

while IFS='' read -r line || [[ -n "$line" ]]; do
    "${SAMTOOLS}" view -F 4 "${BAMFILE}" $line
done < "$REGIONSFILE" 
