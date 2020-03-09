#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=5:mem=6gb

#### remember job numbers 0 to 23 ###
# qsub -J 0-23 map_align.sh
	# the 48 reads are paired
	# python index at 0

# completed
# 0-5 completed 

#----- load modules ----#
echo '=================================='
echo -e "\nPython script\n"

module load anaconda3/personal

echo '=================================='
echo -e "\nLoading BWA\n"
module load bwa/0.7.8


#echo '=================================='
#echo -e "\nCopying script\n"
#cp /rds/general/user/jb1919/home/genomics/code/map_align.py $TMPDIR




REF=$EPHEMERAL/test_align/ref_genome_indexed/EquCab2.fna
read1=$EPHEMERAL/test_align/V300044309_L2_B5GHORlfyRAAAAAAA-517_1.fq
read2=$EPHEMERAL/test_align/V300044309_L2_B5GHORlfyRAAAAAAA-517_2.fq

echo '=================================='
echo -e "\nAlign sequence 1 \n"
bwa aln $REF $read1 > $EPHEMERAL/test_align/read1.sai

echo '=================================='
echo -e "\nAlign sequence 2 \n"

bwa aln $REF $read2 > $EPHEMERAL/test_align/read2.sai

echo '=================================='
echo -e "\nbwa sampe \n"

bwa sampe $REF $EPHEMERAL/test_align/read1.sai $EPHEMERAL/test_align/read2.sai $read1 $read2 > $EPHEMERAL/test_align/aln-pe.sam


