#! /bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=32:mem=62gb

# qsub ancestry.sh

CODE_DIR=$HOME/genomics/code/

cp $CODE_DIR/bamAncestry.sh $CODE_DIR/unix_functions.sh $TMPDIR


echo "Checking Quality"

#bash bamAncestry.sh -q


#echo '=================================='
#echo -e "\nDetermining Chromosome number\n"
#
#if (($PBS_ARRAY_INDEX == 32)); then
#	CHR=chrX
#else
#	CHR=chr$PBS_ARRAY_INDEX	
#fi
#
#
#echo "Chromosome: " $CHR
#echo "Job Number: " $PBS_ARRAY_INDEX

#bash bamAncestry.sh -q chr3
#bash bamAncestry.sh -g $CHR
#bash bamAncestry.sh -c chr3

#bash bamAncestry.sh -p # requires quite a bit of memory
#timer


bash bamAncestry.sh -g 
timer
bash bamAncestry.sh -p
timer