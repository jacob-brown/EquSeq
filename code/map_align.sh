#!/bin/bash
# align pair-ended reads, convert to bam, index, summary stats


DIR=$EPHEMERAL/mapping/


# catch input files 
read REF_GEN BASE_NAME

# create names and paths
FILE_1=$EPHEMERAL/mapping/trimmed/$BASE_NAME
FILE_2

echo '=================================='
echo -e "\nAlign sequences\n"

bwa mem $REF_GEN $FILE_1 $FILE_2 > $DIR/read_out.sam

echo '=================================='
echo -e "\nConvert to bam\n"


samtools view -Sb $DIR/read_out.sam > $DIR/read_out.bam


echo '=================================='
echo -e "\nSorting\n"

samtools sort -m 60GiB  $DIR/read_out.bam -o  $DIR/read_out.sorted.bam

echo '=================================='
echo -e "\nIndex\n"

samtools index $DIR/read_out.sorted.bam

echo '=================================='
echo -e "\nFlagstat\n"

samtools flagstat $DIR/read_out.sorted.bam > $DIR/read_stat.txt
