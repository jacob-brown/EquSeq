#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-04-01
# Last Modified: 2020-05-25
# Desc: Update the bam files with the chromosome location
	# rather than ncbis assecison codes


# sh updateHeader.sh EquCab2.fna 

SCRIPT_DIR=$HOME/genomics/EquSeq/scripts/
# import unix functions
source $SCRIPT_DIR/unix_functions.sh

FILE_IN=$1 # which file to manipulate
REF_TAB=assembly_report.txt 
HEADER=header.bam


echo '=================================='
read -r -p "Is this a reference genome? [y/N] " response


echo '=================================='
echo -e "\nLoading modules\n"
module load samtools/1.3.1 
module load anaconda3/personal


echo '=================================='
echo -e "\nGet assembly report\n"

#wget -O $REF_TAB ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Equus_caballus/all_assembly_versions/suppressed/GCF_000002305.2_EquCab2.0/GCF_000002305.2_EquCab2.0_assembly_report.txt

wget -O $REF_TAB ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Equus_caballus/all_assembly_versions/GCF_002863925.1_EquCab3.0/GCF_002863925.1_EquCab3.0_assembly_report.txt

timer

##############################################
################## RESPONSE ##################
##############################################

response=${response,,}    # to lower
if [[ "$response" =~ ^(yes|y)$ ]]
then

#---------------- fasta ----------------#

echo '=================================='
echo -e "\nUpdate fasta header\n"

python3 $SCRIPT_DIR/updateHeader.py 'fasta' $REF_TAB $FILE_IN

timer

echo '=================================='
echo -e "\nClean up the environment\n"

rm $REF_TAB

timer

else


#---------------- bam ----------------#

echo '=================================='
echo -e "\nCreate bam header\n"
samtools view -H $FILE_IN > $HEADER

timer

echo '=================================='
echo -e "\nUpdating bam header\n"

# update and reheader file
python3 $SCRIPT_DIR/updateHeader.py 'bam' $REF_TAB $HEADER | \
	samtools reheader - $FILE_IN > $FILE_IN.newhead.bam

timer

echo '=================================='
echo -e "\nClean up the environment\n"

rm $REF_TAB $HEADER

timer

fi

################ end RESPONSE ################

echo '=================================='
echo -e "\nFinished. Exiting.\n"



