#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=1:mem=1gb

# path to files
PATH_IN='/rds/general/user/jb1919/home/genomics/sequences/horse/cdts-hk.genomics.cn/Clean/F19FTSEUHT1854-swab-horse-1A/'

# path to output
PATH_OUT='/rds/general/user/jb1919/home/genomics/results/report/'

# list of files
FILES=$(ls $PATH_IN/*.fq.gz | tail -1)

# load fastqc
echo "load fastqc."
module load fastqc/0.11.5

# loop through files and generate report
echo "Generate reports."

for f in $FILES;
    do
    	fastqc -d . -o $PATH_OUT $f 
    done

# end
echo 'Report generation complete'
exit
