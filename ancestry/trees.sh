#! /bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=5:mem=5gb


# qsub -J 0-6 trees.sh

CODE_DIR=~/genomics/EquSeq/
DIR=$EPHEMERAL/ancestry/treemix/
FILE=$DIR/treemix.frq.gz
M=$PBS_ARRAY_INDEX
OUT=$DIR/results/no_jackknife/tree.out.$M
TREEMIX=$EPHEMERAL/dependencies/treemix-1.13/src/treemix 
source $CODE_DIR/scripts/unix_functions.sh

echo '=================================='
echo -e "\nrunning treemix\n"
$TREEMIX -i $FILE -m $M -o $OUT -noss -n_warn 5 -root Przewalski > $OUT.log

timer
# 13m for 0 migrations

# &> /dev/null 


# ~40mins for 5 migrations 


