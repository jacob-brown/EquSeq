#!/bin/bash
# trim fastq files

# catch input files 
READ1=$1
READ2=$2
BASE_NAME=$3

BASE1=$(basename $READ1 | cut -f 1 -d '.')
BASE2=$(basename $READ2 | cut -f 1 -d '.')

FILE_OUT1=$EPHEMERAL/mapping/trimmed/$BASE_NAME'_1.trim.fq'
FILE_OUT2=$EPHEMERAL/mapping/trimmed/$BASE_NAME'_2.trim.fq'


echo -e 'trimming: ' $READ1 'and' $READ2 

# trim
echo '-----------------------'
echo -e "\nTrimming\n"

echo 'read 1'
fastx_trimmer -l 90 -i $READ1 -o $FILE_OUT1
		# ~ 8 mins per read

echo 'read 2'
fastx_trimmer -l 90 -i $READ2 -o $FILE_OUT2

exit


