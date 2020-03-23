#!/bin/bash
# align pair-ended reads, convert to bam, index, summary stats


# catch input files 
BASE_NAME=$1 

DIR=$EPHEMERAL/mapping/
REF_GEN=$DIR/ref_genome/EquCab2.fna

# create names and paths
FILE_1=$DIR/trimmed/$BASE_NAME'_1.trim.fq'
FILE_2=$DIR/trimmed/$BASE_NAME'_2.trim.fq'

echo '-----------------------'
echo -e "\nAlign sequences\n"

bwa mem $REF_GEN $FILE_1 $FILE_2 > $DIR/aligned/$BASE_NAME'.sam'

echo '-----------------------'
echo -e "\nConvert to bam\n"


samtools view -Sb $DIR/aligned/$BASE_NAME'.sam' > $DIR/converted/$BASE_NAME'.bam'


echo '-----------------------'
echo -e "\nSorting\n"

samtools sort -m 60GiB  $DIR/converted/$BASE_NAME'.bam' -o  \
		$DIR/sorted/$BASE_NAME'.sorted.bam'

echo '-----------------------'
echo -e "\nIndex\n"

samtools index $DIR/sorted/$BASE_NAME'.sorted.bam'

echo '-----------------------'
echo -e "\nFlagstat\n"

samtools flagstat $DIR/sorted/$BASE_NAME'.sorted.bam' > \
		$DIR/stats/$BASE_NAME'.stat.txt'
