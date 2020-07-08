#! /bin/bash
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=1:mem=1gb

CODE_DIR=~/genomics/EquSeq/
DIR=$EPHEMERAL/ancestry/treemix/
FILE=$DIR/treemix.frq.gz
OUT=$DIR/results/f3stat.txt
THREEPOP=$EPHEMERAL/dependencies/treemix-1.13/src/threepop 
source $CODE_DIR/scripts/unix_functions.sh

echo '=================================='
echo "\nrun threepop (treemix f3 stat)\n"
# blocks of 10 snps, as in sheep paper 
$THREEPOP -i $FILE -k 500 > $OUT

timer
