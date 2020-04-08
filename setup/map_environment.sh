#!/bin/bash
#PBS -l walltime=03:00:00
#PBS -l select=1:ncpus=1:mem=1gb

# Desc: # Establish the environment for running the cleaning and mapping 

DIR=$EPHEMERAL/mapping/

echo '=================================='
echo -e "\nCreate directories\n"

# create the folder directories
mkdir \
$DIR/aligned/ \
$DIR/converted/ \
$DIR/merged/ \
$DIR/reads/ \
$DIR/ref_genome/ \
$DIR/sorted/ \
$DIR/stats/ \
$DIR/trimmed/ \
$DIR/tmp_cleaning/\
$DIR/cleaned/





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


