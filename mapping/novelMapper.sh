#!/bin/bash
#PBS -l walltime=40:00:00
#PBS -lselect=1:ncpus=30:mem=360gb


# qsub -J 0-23 master.sh


CODE_DIR=$HOME/genomics/EquSeq/
DIR=$EPHEMERAL/novel_data/

# import unix functions
source $HOME/genomics/EquSeq/scripts/unix_functions.sh


echo '=================================='
echo -e "\nLoading modules\n"
module load fastx/0.0.14 # trimming 
module load bwa/0.7.8 # alignment
module load samtools/1.3.1 # general
module load java/jdk-8u144 # picard associated
module load picard/2.6.0 # cleaning
module load anaconda3/personal # python
module load htslib


# file names based on job number 
FILE=($(python3 $CODE_DIR/mapping/names.py | tr -d "[''],"))

# pair ended reads
READ1=${FILE[0]} # 1st pair ended read
READ2=${FILE[1]} # 2nd pair ended read

# new shorter file name
FILE_PREFIX=${FILE[2]}

echo "read1: "$READ1
echo "read2: "$READ2
echo "prefix: "$FILE_PREFIX


echo '=================================='
echo -e "\nCopy reads from HOME\n"

# copy read files to ephemeral
# read 1
cp $HOME/genomics/sequences/cdts-hk.genomics.cn/Clean/F19FTSEUHT1854-swab-horse-1A/$READ1 \
	$DIR/raw_files/$FILE_PREFIX'_1.fq.gz'

# read 2
cp $HOME/genomics/sequences/cdts-hk.genomics.cn/Clean/F19FTSEUHT1854-swab-horse-1A/$READ2 \
	$DIR/raw_files/$FILE_PREFIX'_2.fq.gz'

timer 
echo '=================================='
echo -e "\nUnzipping\n"

gunzip $DIR/raw_files/$FILE_PREFIX'_1.fq.gz'
gunzip $DIR/raw_files/$FILE_PREFIX'_2.fq.gz'

timer 
echo '=================================='
echo -e "\nTrimming\n"

FILE_1=$DIR/trimmed/$BASE_NAME'_1.trim.fq'
FILE_2=$DIR/trimmed/$BASE_NAME'_2.trim.fq'

echo 'read 1'
fastx_trimmer -l 90 -i $DIR/raw_files/$FILE_PREFIX'_1.fq' -o $FILE_1
		# ~ 8 mins per read
echo 'read 2'
fastx_trimmer -l 90 -i $DIR/raw_files/$FILE_PREFIX'_2.fq' -o $FILE_2

timer


echo '=================================='
echo -e "\nAligning\n"


sh $CODE_DIR/mapping/align.sh $FILE_1 $FILE_2 $FILE_PREFIX $DIR

# timer
timer 


