#! /bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=5:mem=5gb


# qsub -J 0-5 trees_bootstrap.sh

CODE_DIR=~/genomics/EquSeq/
DIR=$EPHEMERAL/ancestry/treemix/
FILE=$DIR/treemix.frq.gz
M=$PBS_ARRAY_INDEX
TREEMIX=$EPHEMERAL/dependencies/treemix-1.13/src/treemix 
source $CODE_DIR/scripts/unix_functions.sh

echo '=================================='
echo -e "\nrunning treemix\n"

echo "20 repeats"
for i in {0..20}
do
	OUT=$DIR/results/jackknife/tree.out.M$M.I$i
	$TREEMIX -i $FILE -m $M -k 500 -o $OUT \
		-noss -n_warn 5 -root Przewalski > $OUT.log
	timer
done

timer
# 13m for 0 migrations

# &> /dev/null 

# ~40mins for 5 migrations 

$TREEMIX -i $FILE -m 1 -k 500 -o test -noss -n_warn 5 -root Przewalski > test.log