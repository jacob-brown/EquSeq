#!/bin/bash
#PBS -l walltime=22:00:00
#PBS -l select=1:ncpus=32:mem=62gb


# Run job for 2 files at once _1 and _2
# after the mapping this will make it easier to handle
# i.e. trimming will handle 2 files at once


# qsub -J 0-23 map_master.sh

# import unix functions
source ../scripts/unix_functions.sh

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
cp $HOME/genomics/code/map_names.py \
   $HOME/genomics/code/map_trim.sh \
   $HOME/genomics/code/map_align.sh \
   $TMPDIR


echo '=================================='
echo -e "\nDefine paths and files\n"

DIR=$EPHEMERAL/mapping/

# file names based on job number 
FILE=($(python3 map_names.py | tr -d "[''],"))

# pair ended reads
READ1=${FILE[0]} # 1st pair ended read
READ2=${FILE[1]} # 2nd pair ended read

# new shorter file name
FILE_PREFIX=${FILE[2]}

echo "read1: "$READ1
echo "read2: "$READ2
echo "prefix: "$FILE_PREFIX

# index reference genome
# trim
# align
# convert to bam
# sort
# index 
# stats
# cleaning
	# fixmate
	# remove duplicates
		# check
	# mapped
	# unmapped
	# merge

# timer
timer 

echo '==================================================='
echo '==================================================='
echo -e '\n---------- Begin Running Scripts ----------\n'
echo '==================================================='
echo '==================================================='

echo '=================================='
echo -e "\nCopy reads from HOME\n"

# copy read files to ephemeral
# read 1
cp $HOME/genomics/sequences/cdts-hk.genomics.cn/Clean/F19FTSEUHT1854-swab-horse-1A/$READ1 \
	$DIR/reads/$FILE_PREFIX'_1.fq.gz'

# read 2
cp $HOME/genomics/sequences/cdts-hk.genomics.cn/Clean/F19FTSEUHT1854-swab-horse-1A/$READ2 \
	$DIR/reads/$FILE_PREFIX'_2.fq.gz'

timer 
echo '=================================='
echo -e "\nUnzipping\n"

gunzip $DIR/reads/$FILE_PREFIX'_1.fq.gz'
gunzip $DIR/reads/$FILE_PREFIX'_2.fq.gz'

timer 
echo '=================================='
echo -e "\nTrimming\n"

# trimming reads
	# files in and file out required
bash map_trim.sh $DIR/reads/$FILE_PREFIX'_1.fq' \
	$DIR/reads/$FILE_PREFIX'_2.fq' $FILE_PREFIX

# timer
timer 

echo '=================================='
echo -e "\nAligning\n"

bash map_align.sh $FILE_PREFIX

# timer
timer 








