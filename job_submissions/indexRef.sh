#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=32:mem=62gb


REF=$EPHEMERAL/mapping/ref_genome/EquCab2.fna

source $HOME/genomics/EquSeq/scripts/unix_functions.sh

echo '=================================='
echo -e "\nLoading modules\n"

module load samtools/1.3.1 # general


echo '-----------------------'
echo -e "\nSamtools Index\n"

samtools faidx $REF

timer