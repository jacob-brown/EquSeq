#! /bin/bash
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=1:mem=1gb

#qsub -J 0-99 window.sh 
#qsub -J 0-399 window.sh 
#qsub -J 0-3 window.sh 

CODE_DIR=$HOME/genomics/EquSeq/

source $CODE_DIR/scripts/unix_functions.sh

module load anaconda3/personal

DIR=$EPHEMERAL/ancestry/
ALL_FILES=($( ls -v $DIR/all_beagles/*beagle.gz ))
FILE=${ALL_FILES[$PBS_ARRAY_INDEX]} 
BASE=($( basename $FILE))

echo "beagle file: " $FILE

### beagle to list ### 

python $CODE_DIR/ancestry/snpHandler.py -c window \
	-i $FILE -o $DIR/all_beagles/$BASE.list.out -d 5


### subset ###
zcat $FILE | head -1 > $FILE.TMP.HEAD 

zcat $FILE | grep -Fw -f $DIR/all_beagles/$BASE.list.out - \
	| cat $FILE.TMP.HEAD - > $DIR/filtered_beagles/$BASE.filter.beagle

gzip $DIR/filtered_beagles/$BASE.filter.beagle
rm $FILE.TMP.HEAD

exit

#python $CODE_DIR/ancestry/snpHandler.py -c subBeagle \
#	-i $FILE \
#	-o $DIR/$BASE.filter \
#	-l $DIR/$BASE.list.out

#gzip $DIR/$BASE.filter.beagle