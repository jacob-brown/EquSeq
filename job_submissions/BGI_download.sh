#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=1gb

echo "script starting"

# wget -m ftp://account:password@link_from_bgi
#wget -m ftp://F19FTSEUHT1854:HORlfyR_1@cdts-hk.genomics.cn/readme.txt -P $HOME/genomics/sequences/horse

wget -m ftp://F19FTSEUHT1854:HORlfyR_1@cdts-hk.genomics.cn/readme.txt -P /rds/general/project/savolainen-archive-2018/live/rawdata/Horses

# or for a specific file 
#wget -m ftp://account:password@link_from_bgi/readme.txt 


echo "data complete."
