#!/bin/bash
#PBS -l walltime=06:00:00
#PBS -l select=1:ncpus=1:mem=5gb

DATA=$(ls $EPHEMERAL/reads/*)


echo '=================================='
echo -e "\nLoading fastx\n"
module load fastx/0.0.14

echo '=================================='
echo -e "\nLoading fastqc\n"
module load fastqc/0.11.5

echo '=================================='
echo '=================================='
echo -e "\n       Begin Loop\n"
echo '=================================='
echo '=================================='

for f in $DATA:
	do
		# unzip
		echo '----------------------'
		echo -e "\nUnzipping\n"

		gunzip $f

		echo '----------------------'
		echo -e "\nPreparing dirs\n"

		# prep file extensions
		NOEXT=$(echo $f | cut -f 1 -d '.') # trim 
		BASE=$(basename "$NOEXT") # basename only
		FILE=$NOEXT.fq # after unzip
		FILE_OUT=$EPHEMERAL/trimmed/$BASE.fq # trim noted
		

		# trim
		echo '----------------------'
		echo -e "\nTrimming\n"

		
		fastx_trimmer -l 90 -i $FILE -o $FILE_OUT
		# ~ 8 mins per read

		# generate report
		echo '----------------------'
		echo -e "\nReporting\n"

		fastqc -d . -o $EPHEMERAL/report $FILE_OUT

	done
