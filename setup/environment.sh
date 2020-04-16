#!/bin/bash
#PBS -l walltime=03:00:00
#PBS -l select=1:ncpus=1:mem=1gb

# Desc: # Establish the environment for running the cleaning and mapping 



echo '=================================='
echo -e "\nCreate directories\n"

# novel sample(s)
mkdir $EPHEMERAL/new_data/ &&\
DIR_SAMPLE=$EPHEMERAL/new_data/

# all available data 
mkdir $EPHEMERAL/all_data/ &&\
DIR_ALL=$EPHEMERAL/all_data/

mkdir \
	$DIR_SAMPLE/aligned/ \
	$DIR_SAMPLE/converted/ \
	$DIR_SAMPLE/merged/ \
	$DIR_SAMPLE/reads/ \
	$DIR_SAMPLE/ref_genome/ \
	$DIR_SAMPLE/sorted/ \
	$DIR_SAMPLE/stats/ \
	$DIR_SAMPLE/trimmed/ \
	$DIR_SAMPLE/cleaned/ 

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
echo -e "\nDownload reference genome\n"

# get reference genome
#-------------- EquCab2 --------------#

# specific file 
rsync --copy-links --times --verbose \
	rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Equus_caballus/all_assembly_versions/suppressed/GCF_000002305.2_EquCab2.0/GCF_000002305.2_EquCab2.0_genomic.fna.gz \
			$DIR/ref_genome/


#-------------- EquCab3 --------------#
# specific file 
#rsync --copy-links --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Equus_caballus/latest_assembly_versions/GCF_002863925.1_EquCab3.0/GCF_002863925.1_EquCab3.0_genomic.fna.gz .



echo '=================================='
echo -e "\nExit\n"

exit


