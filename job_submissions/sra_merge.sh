#!/bin/bash
#PBS -l walltime=22:00:00
#PBS -l select=1:ncpus=32:mem=62gb

# qsub -J 0-35 wgs_merge.py.sh

# import unix functions
cp $HOME/genomics/code/unix_functions.sh \
	$HOME/genomics/code/wgs_merge.py.py \
	$TMPDIR

source unix_functions.sh

#----- load modules ----#
echo '=================================='
echo -e "\nLoad samtools\n"
module load samtools/1.3.1
module load anaconda3/personal # python

#----- merge ----#
echo '=================================='
echo -e "\nMerge Files\n"

python wgs_merge.py.py

timer
