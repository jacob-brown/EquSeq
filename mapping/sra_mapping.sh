#!/bin/bash
#PBS -l walltime=46:00:00
#PBS -lselect=1:ncpus=30:mem=360gb

# ERR868003, ERR868004 - memory maxed when aligning

# qsub -J 0-47 sra_mapping.sh
	# complete for sta_runs 
	# 0-1 for sra_runs_to_do

# qsub -J 0-49 sra_mapping.sh
# qsub -J 0-581 sra_mapping.sh


echo '=================================='
echo -e "\nLoading modules\n"
module load fastx/0.0.14 # trimming 
module load bwa/0.7.8 # alignment
module load samtools/1.3.1 # general
module load java/jdk-8u144 # picard associated
module load picard/2.6.0 # cleaning
module load anaconda3/personal # python


echo '=================================='
echo -e "\nMove scripts to TMPDIR\n"
cp $HOME/genomics/code/sra_align.sh \
	$HOME/genomics/code/sra_mapping.py \
	$HOME/genomics/code/unix_functions.sh \
	$TMPDIR

source unix_functions.sh

echo '=================================='
echo -e "\nDefine paths and files\n"

DIR=$EPHEMERAL/all_data/files_to_align/files/
LIST_FASTQ=$EPHEMERAL/all_data/sra_runs_to_do.txt # sra_runs.txt

# file names based on job number 
FILE=($(python3 sra_mapping.py $LIST_FASTQ $DIR | tr -d "[''],"))

# pair ended reads
READ1=${FILE[0]} # 1st pair ended read
READ2=${FILE[1]} # 2nd pair ended read

# new shorter file name
FILE_PREFIX=${FILE[2]}

echo "read1: "$READ1
echo "read2: "$READ2
echo "prefix: "$FILE_PREFIX

timer 

echo '==================================================='
echo '==================================================='
echo -e '\n---------- Begin Running Scripts ----------\n'
echo '==================================================='
echo '==================================================='


echo '=================================='
echo -e "\nAligning\n"

bash sra_align.sh $READ1 $READ2 $FILE_PREFIX

# timer
timer 



