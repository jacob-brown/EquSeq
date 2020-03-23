#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=1:mem=1gb

# Desc: # Establish the environment for running the cleaning and mapping 

DIR=$EPHEMERAL/mapping/
DATA=($DIR/reads/*.gz) # array of all data _1 and _2 extensions
FILE=${DATA[$PBS_ARRAY_INDEX]} # select the data by the job number

echo '=================================='
echo -e "\nUnzip reads\n"

echo -e "\nreads\n"
# reads
gunzip $FILE

echo '=================================='
echo -e "\nExit\n"

exit


