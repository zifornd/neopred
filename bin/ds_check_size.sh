#!/bin/bash
#echo "Printing first parameter to the script: $1"
#cat "$1"
size=$(cat "$1" | awk 'FNR == 8 {print}' | grep -o '[[:digit:]]*')
#echo "Printing size variable inside script: $size"
size=$(($size/1000000))
#echo "Printing size variable 2nd time inside script: $size"
if [ $size -lt 10 ]
then
    downsampling_size=$size'M'
    #s="$size"
else
    downsampling_size=$(($size*2/10))'M'
    #s="$size"
fi
#echo $downsampling_size > $2
#echo "$size"
echo $downsampling_size
