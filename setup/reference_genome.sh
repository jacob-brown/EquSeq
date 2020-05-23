#!/bin/bash
#PBS -l walltime=30:00:00
#PBS -l select=1:ncpus=32:mem=124gb

# Desc: # retrieve and index reference genome

CODE_DIR=$HOME/genomics/EquSeq/

# import unix functions
source $CODE_DIR/scripts/unix_functions.sh

DIR_REF=$EPHEMERAL/ref_genome/


echo '=================================='
echo -e "\nLoading Modules\n"
module load bwa/0.7.8



echo '=================================='
echo -e "\nDownload reference genome\n"

# get reference genome
#-------------- EquCab2 --------------#

# specific file 
#rsync --copy-links --times --verbose \
#	rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Equus_caballus/all_assembly_versions/suppressed/GCF_000002305.2_EquCab2.0/GCF_000002305.2_EquCab2.0_genomic.fna.gz \
			#$DIR_REF


#-------------- EquCab3 --------------#
# specific file 
rsync --copy-links --times --verbose \
	rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Equus_caballus/latest_assembly_versions/GCF_002863925.1_EquCab3.0/GCF_002863925.1_EquCab3.0_genomic.fna.gz $DIR_REF


timer

echo '=================================='
echo -e "\nUnzipping reference genome\n"

gunzip $DIR_REF/*.gz 

echo '=================================='
echo -e "\nRename\n"

#mv $DIR_REF/GCF_000002305.2_EquCab2.0_genomic.fna $DIR_REF/EquCab2.fna
mv  $DIR_REF/GCF_002863925.1_EquCab3.0_genomic.fna $DIR_REF/EquCab3.fna

echo '=================================='
echo -e "\nIndex ref genome\n"


#bwa index $DIR_REF/EquCab2.fna
bwa index $DIR_REF/EquCab3.fna
# roughly 30 mins


echo '=================================='
echo -e "\nExit\n"

exit



