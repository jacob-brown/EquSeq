#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=10gb

echo "script starting"

# wget -m ftp://account:password@link_from_bgi
wget -m ftp://F19FTSEUHT1854:HORlfyR_1@cdts-hk.genomics.cn/readme.txt -P $HOME/genomics/sequences/horse

# or for a specific file 
#wget -m ftp://account:password@link_from_bgi/readme.txt 

echo "download complete"

# move all the data to the desied folder in home dir
mv * $HOME/HorseSequence

echo "data moved."
