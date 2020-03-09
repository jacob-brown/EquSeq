#!/bin/bash
#PBS -l walltime=02:00:00
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


FILE=$EPHEMERAL/test_align/aln-pe.sam

#----- converting to sam ----#
FILE_CONV_OUT=$EPHEMERAL/test_align/aln-pe_con.bam # outfile of bam to sam

#----- sorting sam ----#
FILE_SORT_OUT=$EPHEMERAL/test_align/aln-pe_con.sorted.bam # outfile of sorted sam


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