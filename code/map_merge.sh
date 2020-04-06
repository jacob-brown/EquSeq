#!/bin/bash
#PBS -l walltime=22:00:00
#PBS -l select=1:ncpus=32:mem=62gb


#DATA=($EPHEMERAL/mapping/cleaned/*.bam) # array of all data
DATA=($EPHEMERAL/mapping/sorted/*.bam)
FILE_OUT=($EPHEMERAL/mapping/merged/merged_reads.bam)

#----- load modules ----#
echo '=================================='
echo -e "\nLoad samtools\n"
module load samtools/1.3.1

#----- merge ----#
echo '=================================='
echo -e "\nMerge Files\n"

samtools merge -r --threads 32 $FILE_OUT $DATA

echo '=================================='
echo -e "\nIndex\n"

samtools index --threads 32 $FILE_OUT $FILE_OUT.bai


echo '=================================='
echo -e "\nflagstat\n"

samtools flagstat $FILE_OUT > \
		$EPHEMERAL/merged/merged_reads'.stat.txt'
