#!/bin/bash
#PBS -l walltime=03:00:00
#PBS -l select=1:ncpus=1:mem=1gb

# Desc: # Establish the environment for running the cleaning and mapping 

module load anaconda3/personal # python

echo '=================================='
echo -e "\nCreate directories\n"

# novel sample(s)
mkdir $EPHEMERAL/new_data/ &&\
DIR_SAM=$EPHEMERAL/new_data/

# reference genome
mkdir $EPHEMERAL/ref_genome/ &&\
DIR_REF=$EPHEMERAL/ref_genome/

# all available data 
mkdir $EPHEMERAL/all_data/ &&\
DIR_ALL=$EPHEMERAL/all_data/

# ancestry data store
mkdir $EPHEMERAL/ancestry/
mkdir $EPHEMERAL/ancestry/pca
mkdir $EPHEMERAL/ancestry/qualityChecks

# gene to trait data store
mkdir $EPHEMERAL/gene_to_trait/

# dependencies
mkdir $EPHEMERAL/dependencies/ &&\
DIR_DEP=$EPHEMERAL/dependencies/

mkdir \
	$DIR_SAM/aligned/ \
	$DIR_SAM/converted/ \
	$DIR_SAM/merged/ \
	$DIR_SAM/reads/ \
	$DIR_SAM/ref_genome/ \
	$DIR_SAM/sorted/ \
	$DIR_SAM/stats/ \
	$DIR_SAM/trimmed/ \
	$DIR_SAM/cleaned/ 

mkdir \
	$DIR_ALL/aligned/ \
	$DIR_ALL/converted/ \
	$DIR_ALL/files/ \
	$DIR_ALL/sorted/ \
	$DIR_ALL/stats/ \
	$DIR_ALL/cleaned/ 
	#$DIR_ALL/merged/ \
	##$DIR_ALL/trimmed/ \
	#$DIR_ALL/ref_genome/ \


echo '=================================='
echo -e "\nBuild angsd\n"


wget http://popgen.dk/software/download/angsd/angsd0.930.tar.gz 
tar xf angsd0.930.tar.gz
cd htslib;make;cd ..
cd angsd
make HTSSRC=../htslib
cd ..
rm angsd0.930.tar.gz

echo '=================================='
echo -e "\nBuild pcangsd\n"


git clone https://github.com/Rosemeis/pcangsd.git $DIR_DEP/pcangsd

python $DIR_DEP/pcangsd/setup.py build_ext --inplace

# check if functional
# python pcangsd.py -h


echo '=================================='
echo -e "\nExit\n"

exit


