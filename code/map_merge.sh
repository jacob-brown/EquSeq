#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=5:mem=5gb


DATA=($EPHEMERAL/sorted/*.bam) # array of all data
FILE_OUT=($EPHEMERAL/merged/merged_reads.bam)

#----- load modules ----#
echo '=================================='
echo -e "\nLoad samtools\n"
module load samtools/1.3.1


#----- merge ----#
echo '=================================='
echo -e "\nMerge Files\n"


samtools merge $FILE_OUT $DATA

echo '=================================='
echo -e "\nIndex\n"

samtools index $FILE_OUT