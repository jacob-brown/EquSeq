#! /bin/bash
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=5:mem=5gb

# qsub pca.sh 

module load anaconda3/personal
CODE_DIR=$HOME/genomics/EquSeq/
PCANGSD=$EPHEMERAL/dependencies/pcangsd/pcangsd.py
DIR=$EPHEMERAL/ancestry/
source $CODE_DIR/scripts/unix_functions.sh

source activate myenv # issues with install

echo '=================================='
echo -e "\nall breeds pca 0.02\n"

python $PCANGSD -beagle $DIR/ALL.merged.beagle.gz \
	-o $DIR/ALL.02.PCA -threads 4 -minMaf 0.02 -sites_save

timer

echo '=================================='
echo -e "\nall breeds pca 0.05\n"
python $PCANGSD -beagle $DIR/ALL.merged.beagle.gz \
	-o $DIR/ALL.05.PCA -threads 4 -minMaf 0.05 -sites_save

timer


echo '=================================='
echo -e "\nno prw pca\n"

# default -minMaf 0.05
python $PCANGSD -beagle $DIR/NO.PRZ.merged.beagle.gz \
	-o $DIR/NO.PRZ.PCA -threads 4 -minMaf 0.02 -sites_save


timer
conda deactivate # issues with install

