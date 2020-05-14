#!/bin/bash
#PBS -l walltime=22:00:00
#PBS -l select=1:ncpus=32:mem=62gb

# import unix functions
source ../scripts/unix_functions.sh

#DATA=($EPHEMERAL/mapping/cleaned/*.bam) # array of all data
DATA=($EPHEMERAL/mapping/sorted/*.bam)
FILE_OUT=$EPHEMERAL/mapping/merged/merged.bam

#----- load modules ----#
echo '=================================='
echo -e "\nLoad samtools\n"
module load samtools/1.3.1

#----- merge ----#
echo '=================================='
echo -e "\nMerge Files\n"

samtools merge --threads 31 $FILE_OUT $DATA

timer

echo '=================================='
echo -e "\nIndex\n"

samtools index $FILE_OUT $FILE_OUT.bai


timer

echo '=================================='
echo -e "\nflagstat\n"

samtools flagstat $FILE_OUT > $EPHEMERAL/mapping/merged/merged'.stat.txt'

