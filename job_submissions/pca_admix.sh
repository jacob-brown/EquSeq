#! /bin/bash
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=1:mem=gb

CODE_DIR=$HOME/genomics/EquSeq/
DIR=$EPHEMERAL/ancestry/

source $CODE_DIR/scripts/unix_functions.sh

module load anaconda3/personal

echo "pca"

# zcat ALL.merged.beagle.gz | wc -l 
# 247611
# zcat ALL.merged.beagle.gz | head -100000 > tmp.beagle && gzip tmp.beagle
# sh $CODE_DIR/ancestry/bamAncestry.sh -p tmp.beagle.gz
# sh $CODE_DIR/ancestry/bamAncestry.sh -a tmp.beagle.gz


# sh $CODE_DIR/ancestry/bamAncestry.sh -p $DIR/ALL.merged.beagle.gz
# 
# timer
# 
# echo "admix"
# 
# sh $CODE_DIR/ancestry/bamAncestry.sh -a $DIR/ALL.merged.beagle.gz
# 
# timer


# qsub -J 1-6 pca_admix.sh 

NGSADMIX=$EPHEMERAL/dependencies/angsd/misc/NGSadmix

K=${ALL_FILES[$PBS_ARRAY_INDEX]} 

$NGSADMIX -likes $DIR/ALL.merged.beagle.gz -K $K\
			-outfiles $DIR/ALL.MIX.K$K -P 4 -minMaf 0.05