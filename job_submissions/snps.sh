#! /bin/bash
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=5:mem=5gb

# uses very little memory, but takes a long time 
	# check for threading

#qsub -J 0-number of files snps.sh 
# qsub -J 0-213 snps.sh 
	# qsub -J 0-3 snps.sh 
	# qsub -J 0-213 snps.sh 

# import unix functions
source $HOME/genomics/EquSeq/scripts/unix_functions.sh
CODE_DIR=$HOME/genomics/EquSeq/
SNP_DIR=$CODE_DIR/data/snp_calling_list/

ALL_FILES=($(ls $SNP_DIR))
FILE=$SNP_DIR${ALL_FILES[$PBS_ARRAY_INDEX]} 
echo $FILE

echo "running snpCaller"

sh $CODE_DIR/scripts/snpCaller.sh $FILE $PBS_ARRAY_INDEX

timer


# old methods
#ALL_CHR=( chr{1..31} )
#ALL_CHR+=(chrM chrX)
#echo "${ALL_CHR[*]}"
#CHR="${ALL_CHR[1]}"
#CHR="${ALL_CHR[$PBS_ARRAY_INDEX]}"
#sh $CODE_DIR/scripts/snpCaller.sh $CHR