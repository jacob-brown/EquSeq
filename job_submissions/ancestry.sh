#! /bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=32:mem=62gb


# gl 
	# qsub -J 0-32 ancestry.sh 
# pca
	# qsub ancestry.sh


#CODE_DIR=$HOME/genomics/code/
CODE_DIR=$HOME/genomics/HorseGenomics/

source $CODE_DIR/scripts/unix_functions.sh

#echo "Checking Quality"

#bash bamAncestry.sh -q


#### 
## GL 
	# break down into chromosomes, requires an array job 
		# of length snp.chr/*
		# qsub -J 0-32 ancestry.sh 
			# wt: 24hr, 1, 32, 124
FILES=($( ls -v $CODE_DIR/data/ancestry/snp.chr/* ))

# pair ended reads $PBS_ARRAY_INDEX
CHRFILE=${FILES[$PBS_ARRAY_INDEX]} # 1st pair ended read

echo $CHRFILE
echo "GL"
sh $CODE_DIR/ancestry/bamAncestry.sh -g $CHRFILE


####
### Merge beagle files 
	# pass directory of in and out files
#sh $CODE_DIR/ancestry/beagleMerged.sh $EPHEMERAL/ancestry/
#
##timer
#echo "PCA"
#sh $CODE_DIR/ancestry/bamAncestry.sh -p $EPHEMERAL/ancestry/ALL.merged.beagle.gz
#timer