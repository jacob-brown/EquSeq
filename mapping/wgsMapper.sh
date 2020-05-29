#!/bin/bash
#PBS -l walltime=46:00:00
#PBS -lselect=1:ncpus=30:mem=360gb

# ERR868003, ERR868004 - memory maxed when aligning

# gb and IE
# qsub -J 0-47 wgsMapper.sh
# rest of the data
# qsub -J 0-581 wgsMapper.sh 


CODE_DIR=$HOME/genomics/EquSeq/
DIR=$EPHEMERAL/wgs_data/
RUNS_LIST=$CODE_DIR/data/cleaned_data/sra_runs.txt
	# sra_runs.txt - gb and ie
	# sra_runs_to_do.txt - all others

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



# file names based on job number 
FILE=($(python3 $CODE_DIR/mapping/wgs_mapping.py $RUNS_LIST $DIR/raw_files/ | tr -d "[''],"))

# pair ended reads
FILE_1=${FILE[0]} # 1st pair ended read
FILE_2=${FILE[1]} # 2nd pair ended read

# new shorter file name
FILE_PREFIX=${FILE[2]}

echo "read1: "$FILE_1
echo "read2: "$FILE_2
echo "prefix: "$FILE_PREFIX

timer 


echo '=================================='
echo -e "\nAligning\n"

FILE_1=$DIR/raw_files/$FILE_1
FILE_2=$DIR/raw_files/$FILE_2

sh $CODE_DIR/mapping/align.sh $FILE_1 $FILE_2 $FILE_PREFIX $DIR
# timer
timer 



