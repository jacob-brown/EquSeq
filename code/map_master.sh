#!/bin/bash
#PBS -l walltime=06:00:00
#PBS -l select=1:ncpus=32:mem=62gb

# Run job for 2 files at once _1 and _2
# after the mapping this will make it easier to handle
# i.e. trimming will handle 2 files at once

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
mv $HOME/genomics/code/map_* $TMPDIR


echo '=================================='
echo -e "\nDefine paths and files\n"

#PICARD=$PICARD_HOME/picard.jar
DIR=$EPHEMERAL/mapping/

# file names based on job number 
FILE=($(python3 map_names.py | tr -d "[''],"))

# pair ended reads
READ1=${FILE[0]} # 1st pair ended read
READ2=${FILE[1]} # 2nd pair ended read

# new shorter file name
FILE_PREFIX=${FILE[2]}

#echo $READ1
#echo $READ2
#echo $FILE_PREFIX

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


echo '==================================================='
echo '==================================================='
echo -e '\n---------- Begin Running Scripts ----------\n'
echo '==================================================='
echo '==================================================='


echo '=================================='
echo -e "\nTrimming\n"

# trimming reads
	# files in and file out required
bash map_trim.sh $DIR/reads/$READ1 $DIR/reads/$READ2 $FILE_PREFIX

echo '=================================='
echo -e "\nAligning\n"

bash map_align.sh $FILE_PREFIX


echo '=================================='
echo -e "\nCleaning\n"

bash map_clean.sh $FILE_PREFIX









