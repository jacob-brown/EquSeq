#! /bin/bash
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=1:mem=1gb

CODE_DIR=~/genomics/EquSeq/
DIR=$EPHEMERAL/ancestry/treemix/
FILE=$DIR/treemix.frq.gz
OUT3=$DIR/results/f3stat.txt
OUT4=$DIR/results/f4stat.txt
THREEPOP=$EPHEMERAL/dependencies/treemix-1.13/src/threepop 
FOURPOP=$EPHEMERAL/dependencies/treemix-1.13/src/fourpop
source $CODE_DIR/scripts/unix_functions.sh

echo '=================================='
echo -e "\nrun threepop (treemix f3 stat)\n"

$THREEPOP -i $FILE -k 500 > $OUT3
timer

echo '=================================='
echo -e "\nrun fourpop (treemix f4 stat)\n"

$FOURPOP -i $FILE -k 500 > $OUT4
timer
