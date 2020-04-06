#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-04-01
# Last Modified: 2020-04-02
# Desc: Update the bam files with the chromosome location
	# rather than ncbis assecison codes

###########################################
# timer function

time_start=$SECONDS

function timer {
	duration=$(($SECONDS - $time_start))
	echo -e "\n..........................\n"
 	echo "Time elapsed: " $duration " sec"
 	echo -e "\n..........................\n"
}


###########################################

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

wget -O $REF_TAB ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Equus_caballus/all_assembly_versions/suppressed/GCF_000002305.2_EquCab2.0/GCF_000002305.2_EquCab2.0_assembly_report.txt

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

python3 correct_header.py 'fasta' $REF_TAB $FILE_IN

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
python3 correct_header.py 'bam' $REF_TAB $HEADER | \
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



