#! /bin/bash
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=8:mem=10gb

# nthread (avg) and memory is 1 or less

# gl 
	# qsub -J 0-399 ancestry.sh # 400 bins 
	# qsub -J 0-99 ancestry.sh  # 100 bins
# bcf 
	#qsub -J 0-49 ancestry.sh  # 50 bins

CODE_DIR=$HOME/genomics/EquSeq/

source $CODE_DIR/scripts/unix_functions.sh


#### 
## GL 

#ALL_FILES=($( ls -v $EPHEMERAL/snp_calling/snp_lists/* ))
#FILE=${ALL_FILES[$PBS_ARRAY_INDEX]} 
#
#echo "snp file: " $FILE
#
#sh $CODE_DIR/ancestry/bamAncestry.sh -g $FILE

#### 
## vcf generation
ALL_FILES=($( ls -v $CODE_DIR/data/ancestry/snp.vcf.list/* ))
FILE=${ALL_FILES[$PBS_ARRAY_INDEX]} 
echo "snp file: " $FILE
sh $CODE_DIR/ancestry/bamAncestry.sh -b $FILE
