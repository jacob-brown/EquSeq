#!/bin/bash
#PBS -l walltime=22:00:00
#PBS -l select=1:ncpus=32:mem=62gb

# qsub -J 0-35 wgs_merge.py.sh

# import unix functions
CODE_DIR=$HOME/genomics/EquSeq/

source $CODE_DIR/scripts/unix_functions.sh


#----- load modules ----#
echo '=================================='
echo -e "\nLoad samtools\n"
module load samtools/1.3.1
module load anaconda3/personal # python

#----- merge ----#
echo '=================================='
echo -e "\nMerge Files\n"

python $CODE_DIR/mapping/wgs_merge.py

timer
