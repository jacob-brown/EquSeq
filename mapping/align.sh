#!/bin/bash
# align pair-ended reads, convert to bam, index, summary stats


REF_GEN=$EPHEMERAL/ref_genome/EquCab3.fna

# create names and paths
FILE_1=$1
FILE_2=$2
BASE_NAME=$3
DIR=$4

echo '-----------------------'
echo -e "\nAligning, converting bam\n"


bwa mem -t 29 $REF_GEN $FILE_1 $FILE_2 | samtools view -bS --threads 29 - > \
		$DIR/converted/$BASE_NAME'.bam'


echo '-----------------------'
echo -e "\nSorting\n"

samtools sort -m 300GiB --threads 29 $DIR/converted/$BASE_NAME'.bam' -o  \
		$DIR/sorted/$BASE_NAME'.sorted.bam'


echo '-----------------------'
echo -e "\nIndex\n"

samtools index $DIR/sorted/$BASE_NAME'.sorted.bam'

echo '-----------------------'
echo -e "\nFlagstat\n"

# should be high due to 
samtools flagstat $DIR/sorted/$BASE_NAME'.sorted.bam' > \
		$DIR/stats/$BASE_NAME'.stat.txt'


