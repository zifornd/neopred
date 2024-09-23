#!/bin/bash
##Script to downsample the reads for RSeQC

echo "Printing first parameter: $1"
echo "Printing second parameter: $2"
VALUE=$(echo $2 | sed -r 's/M//g')
#echo "Printing VALUE variable: $VALUE"
FACTOR=$(samtools idxstats $1 | cut -f3 | awk -v COUNT=$((VALUE*1000000)) 'BEGIN {total=0} {total += $1} END {print COUNT/total}')


echo "Printing VALUE variable: $VALUE"
echo "Printing FACTOR variable: $FACTOR"
if [[ $FACTOR > 1 ]]
    then
    echo '[ERROR]: Requested number of reads exceeds total read count in' $1 '-- exiting' && exit 1
fi

samtools view -s $FACTOR -b $1 > $1"_"$2"_downsampling.bam"
