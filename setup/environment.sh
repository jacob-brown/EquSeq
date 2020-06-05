#!/bin/bash
#PBS -l walltime=03:00:00
#PBS -l select=1:ncpus=1:mem=1gb

# Desc: # Establish the environment for running the cleaning and mapping 

module load anaconda3/personal # python

echo '=================================='
echo -e "\nCreate directories\n"

# novel sample(s)
mkdir $EPHEMERAL/novel_data/ &&\
DIR_NOV=$EPHEMERAL/novel_data/

# reference genome
mkdir $EPHEMERAL/ref_genome/ &&\
DIR_REF=$EPHEMERAL/ref_genome/

# all available data 
mkdir $EPHEMERAL/wgs_data/ &&\
WGS_DATA=$EPHEMERAL/wgs_data/

# ancestry data store
mkdir $EPHEMERAL/ancestry/
mkdir $EPHEMERAL/ancestry/qualityChecks
mkdir $EPHEMERAL/ancestry/filtered_beagles
mkdir $EPHEMERAL/ancestry/all_beagles

# gene to trait data store
mkdir $EPHEMERAL/gene_to_trait/

# gene to snps
mkdir $EPHEMERAL/snp_calling/
mkdir $EPHEMERAL/snp_calling/snp_lists/

#
mkdir $EPHEMERAL/oral_diversity

# dependencies
mkdir $EPHEMERAL/dependencies/ &&\
DIR_DEP=$EPHEMERAL/dependencies/

# snp calling
mkdir $EPHEMERAL/snp_calling

mkdir \
	$DIR_NOV/raw_files/ \
	$DIR_NOV/trimmed/ \
	$DIR_NOV/converted/ \
	$DIR_NOV/sorted/ \
	$DIR_NOV/stats/ \
	$DIR_NOV/merged/ 

mkdir \
	$WGS_DATA/raw_files/ \
	$WGS_DATA/converted/ \
	$WGS_DATA/sorted/ \
	$WGS_DATA/stats/ \
	$WGS_DATA/merged/




echo '=================================='
echo -e "\nExit\n"

exit


