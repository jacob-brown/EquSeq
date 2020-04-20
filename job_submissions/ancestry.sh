#! /bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=32:mem=62gb

CODE_DIR=$HOME/genomics/code/

cp $CODE_DIR/bamAncestry.sh $CODE_DIR/unix_functions.sh $TMPDIR


echo "Checking Quality"

bash bamAncestry.sh -q