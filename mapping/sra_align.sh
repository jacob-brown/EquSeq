#!/bin/bash
# align pair-ended reads, convert to bam, index, summary stats

# import unix functions
source unix_functions.sh


DIR=$EPHEMERAL/all_data/files_to_align/
REF_GEN=$EPHEMERAL/mapping/ref_genome/EquCab2.fna

# create names and paths
FILE_1=$DIR/files/$1
FILE_2=$DIR/files/$2
BASE_NAME=$3

#echo '-----------------------'
#echo -e "\nAlign sequences\n"

#bwa mem $REF_GEN $FILE_1 $FILE_2 -t 7 > $DIR/aligned/$BASE_NAME'.sam'

#timer
echo '-----------------------'
echo -e "\nAligning, sorting, converting bam\n"

bwa mem -t 7 $REF_GEN $FILE_1 $FILE_2 | samtools sort -m 1G -@ 4 -O BAM -o $DIR/sorted/$BASE_NAME'.sorted.bam' -
#bwa mem $REF_GEN $FILE_1 $FILE_2 | samtools view -bS - > out.bam
timer

#echo '-----------------------'
#echo -e "\nConvert to bam\n"


#samtools view -bS --threads 29 $DIR/aligned/$BASE_NAME'.sam' > \
#		$DIR/converted/$BASE_NAME'.bam'
#
#timer
#
#echo '-----------------------'
#echo -e "\nSorting\n"
#
#samtools sort -m 60GiB --threads 29 $DIR/converted/$BASE_NAME'.bam' -o  \
#		$DIR/sorted/$BASE_NAME'.sorted.bam'
#
#timer
#
#echo '-----------------------'
#echo -e "\nIndex\n"
#
#samtools index $DIR/sorted/$BASE_NAME'.sorted.bam'
#
#timer
#
#echo '-----------------------'
#echo -e "\nFlagstat\n"
#
## should be high due to 
#samtools flagstat $DIR/sorted/$BASE_NAME'.sorted.bam' > \
#		$DIR/stats/$BASE_NAME'.stat.txt'
#
#timer
#