#!/bin/bash
#PBS -l walltime=10:00:00
#PBS -l select=1:ncpus=5:mem=6gb

# jobs from 0-N
	# bash index starts at 0
	# 0-5 complete 
		# 6-23


#----- load modules ----#
echo '=================================='
echo -e "\nLoad samtools\n"
module load samtools/1.3.1

#----- Prep Data ----#
echo '=================================='
echo -e "\nPreparing the data\n"
# data

DATA=($EPHEMERAL/aligned/*) # array of all data
FILE=${DATA[$PBS_ARRAY_INDEX]} # select the data by the job number
NOEXT=$(echo $FILE | cut -f 1 -d '.') # remove extension 
BASE=$(basename "$NOEXT") # remove path

#----- converting to sam ----#
FILE_CONV_OUT=$EPHEMERAL/converted/$BASE.bam # outfile of bam to sam

#----- sorting sam ----#
FILE_SORT_OUT=$EPHEMERAL/sorted/$BASE.sorted.bam # outfile of sorted sam


echo '=================================='
echo -e "\nConvert to sam\n"

#samtools view -S -b sample.sam > sample.bam
samtools view -S -b $FILE > $FILE_CONV_OUT


echo '=================================='
echo -e "\nSorting\n"


samtools sort -m 5GiB $FILE_CONV_OUT -o $FILE_SORT_OUT

echo '=================================='
echo -e "\nIndex\n"

samtools index $FILE_SORT_OUT


# then merge
# check with  samtools view -H
# samtools flagstat
	# samtools flagstat myFile.bam > myFile.stats 

# check the times of last run, delete and run a trial