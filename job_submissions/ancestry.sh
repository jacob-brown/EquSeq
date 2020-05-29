#! /bin/bash
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=8:mem=10gb

# nthread (avg) and memory is 1 or less

# gl 
	# qsub -J 0-99 ancestry.sh 
		#qsub -J 0-3 ancestry.sh
		# qsub -J 4-99 ancestry.sh

CODE_DIR=$HOME/genomics/EquSeq/

source $CODE_DIR/scripts/unix_functions.sh


#### 
## GL 

ALL_FILES=($( ls -v $EPHEMERAL/snp_calling/snp_lists/* ))
FILE=${ALL_FILES[$PBS_ARRAY_INDEX]} 

echo "snp file: " $FILE


sh $CODE_DIR/ancestry/bamAncestry.sh -g $FILE


