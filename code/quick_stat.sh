#!/bin/bash
#PBS -l walltime=06:00:00
#PBS -l select=1:ncpus=1:mem=1gb

module load samtools/1.3.1 
DIR=$EPHEMERAL/mapping
DATA=($DIR/converted/*)

for i in "${DATA[@]}"
do
	echo $i ;
	BASE=$(basename $i);
	samtools flagstat $i > $DIR/stats/$BASE'.stat.txt';
done