#!/bin/bash
#PBS -lwalltime=48:00:00
#PBS -lselect=1:ncpus=1:mem=1gb

# Desc: Get files from ncbi sra

###########################################
# timer function

time_start=$SECONDS

function timer {
	duration=$(($SECONDS - $time_start))
	echo -e "\n..........................\n"
 	echo "Time elapsed: " $duration " sec"
 	echo -e "\n..........................\n"
}

###########################################

DIR=$EPHEMERAL/sra_data/


#----- load modules ----#
echo '=================================='
echo -e "\nLoad modules\n"
module load sra-toolkit/2.8.1
module load samtools/1.3.1 

echo '=================================='
echo -e "\nClear the environment\n"

rm -f $HOME/ncbi/public/sra/* # force for no errors

echo '=================================='
echo -e "\nFetching\n"

# run code
#prefetch ERR868003
sam-dump ERR868003 | samtools view -bS - > $DIR/ERR868003.bam

timer

ls $HOME/ncbi/public/sra/*
ls $DIR

#echo '=================================='
#echo -e "\nMoving\n"
#mv * $EPHEMERAL/sra_data/
#prefetch --type bam ERR868004 
#prefetch ERR868004

#fastq-dump -X 5 -Z ERR868003

#prefetch --option-file SraAccList.txt
