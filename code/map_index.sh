#!/bin/bash
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=1:mem=5gb



#----- copy to temp ----#
echo '=================================='
echo -e "\nmoving ref genome to temp\n"

# ref genome
rsync -rltv $HOME/genomics/sequences/refgenome/EquCab2.fna $TMPDIR/EquCab2.fna 


#----- Load Modules ----#
# load BWA 
echo '=================================='
echo -e "\nLoading BWA\n"
module load bwa/0.7.8

#----- index ----#
# index reference genome
echo '=================================='
echo -e "\nIndex ref genome\n"
bwa index EquCab2.fna
# roughly 30 mins

# move indexed ref genome and rename
echo '=================================='
echo -e "\nmoving indexed ref gen\n"
mv EquCab2.* $HOME/genomics/sequences/refgenome/indexed/





