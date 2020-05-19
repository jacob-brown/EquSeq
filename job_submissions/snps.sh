#! /bin/bash
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=32:mem=124gb


#qsub -J 0-32 snps.sh 
#qsub -J 0-2 snps.sh 
#qsub -J 3-32 snps.sh 
# import unix functions
source $HOME/genomics/EquSeq/scripts/unix_functions.sh
CODE_DIR=$HOME/genomics/EquSeq/


ALL_CHR=( chr{1..31} )
ALL_CHR+=(chrM chrX)
#echo "${ALL_CHR[*]}"
#CHR="${ALL_CHR[1]}"
CHR="${ALL_CHR[$PBS_ARRAY_INDEX]}"


echo "running snpCaller"

sh $CODE_DIR/ancestry/snpCaller.sh $CHR

