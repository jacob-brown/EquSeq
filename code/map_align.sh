#!/bin/bash
#PBS -l walltime=06:00:00
#PBS -l select=1:ncpus=1:mem=5gb

#### remember job numbers 0 to 23 ###
# qsub -J 0-23 map_align.sh
	# the 48 reads are paired
	# python index at 0


#----- load modules ----#
echo '=================================='
echo -e "\nPython script\n"

module load anaconda3/personal

echo '=================================='
echo -e "\nLoading BWA\n"
module load bwa/0.7.8


echo '=================================='
echo -e "\nCopying script\n"
cp /rds/general/user/jb1919/home/genomics/code/map_align.py $TMPDIR


echo '=================================='
echo -e "\nAlign sequences\n"

python3 map_align.py 

