#!/bin/bash
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=1:mem=5gb

# index the reference genome

# ref genome
DIR_REF=$EPHEMERAL/mapping/ref_genome/

# load
echo '=================================='
echo -e "\nLoading BWA\n"
module load bwa/0.7.8


echo '=================================='
echo -e "\nUnzipping reference genome\n"

#gunzip $DIR_REF/GCF_000002305.2_EquCab2.0_genomic.fna.gz > \
#		$DIR_REF/EquCab2.fna

# index reference genome
echo '=================================='
echo -e "\nIndex ref genome\n"
bwa index $DIR_REF/EquCab2.fna
# roughly 30 mins



echo '=================================='
echo -e "\nExit\n"

exit




