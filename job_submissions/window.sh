#! /bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=1gb

CODE_DIR=$HOME/genomics/EquSeq/

source $CODE_DIR/scripts/unix_functions.sh

module load anaconda3/personal

DIR=$EPHEMERAL/ancestry/
ALL_FILES=($( ls -v $DIR/*beagle.gz ))
FILE=${ALL_FILES[$PBS_ARRAY_INDEX]} 
BASE=($( basename $FILE))

echo "beagle file: " $FILE

### beagle to list ### 

python3 $CODE_DIR/ancestry/snpHandler.py -c window \
	-i $FILE -o $DIR/$BASE.list -d 5


python3 $CODE_DIR/ancestry/snpHandler.py -c subBeagle \
	-i $FILE \
	-o $DIR/$BASE.filter \
	-l $DIR/$BASE.list