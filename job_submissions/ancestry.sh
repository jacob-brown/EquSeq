#! /bin/bash
#PBS -l walltime=30:00:00
#PBS -l select=1:ncpus=32:mem=62gb


# gl 
	# qsub -J 0-31 ancestry.sh 
# pca
	# qsub ancestry.sh


#CODE_DIR=$HOME/genomics/code/
CODE_DIR=$HOME/genomics/EquSeq/

source $CODE_DIR/scripts/unix_functions.sh


#### 
## GL 
	# break down into chromosomes, requires an array job 
		# of length snp.chr/*
		# qsub -J 0-32 ancestry.sh 
			# wt: 24hr, 1, 32, 124

#FILES=($( ls -v $CODE_DIR/data/ancestry/snp.chr/* ))
FILES=($( ls -v $EPHEMERAL/snp_calling/snp_lists/* ))
CHRFILE=${FILES[$PBS_ARRAY_INDEX]} 

echo $CHRFILE
echo "GL"
sh $CODE_DIR/ancestry/bamAncestry.sh -g $CHRFILE


