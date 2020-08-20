#! /bin/bash
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=1:mem=1gb



### Local

time samtools coverage data/processed_sequences/benson/final.bam > results/coverage/coverage_benson.csv


# qsub -J 0-170 coverage.sh

### HPC
	#conda search samtools
	#conda create -n samtools1.10 samtools=1.10
	#cd ../projects/savolainen-archive-2018/live/rawdata/Horses/bam_files/


# setup


module load anaconda3/personal
source activate samtools1.10

BAMDIR=/rds/general/user/jb1919/projects/savolainen-archive-2018/live/rawdata/Horses/bam_files
ALL_FILES=($( ls -v $BAMDIR/*.bam ))
FILE=${ALL_FILES[$PBS_ARRAY_INDEX]} 
BASE=($( basename $FILE))
echo "file: " $FILE
echo "base: " $BASE


samtools coverage $FILE > $HOME/genomics/stats/$BASE.csv

conda deactivate

