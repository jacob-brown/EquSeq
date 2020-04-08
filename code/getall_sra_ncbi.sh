#!/bin/bash
#PBS -lwalltime=24:00:00
#PBS -lselect=1:ncpus=1:mem=10gb

# Desc: Get files from ncbi sra

# qsub -J 0-47 getall_sra_ncbi.sh

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
RES_DIR=$DIR/files/


#----- load modules ----#
echo '=================================='
echo -e "\nLoad modules\n"
#module load sra-toolkit/2.8.1
module load anaconda3/personal



echo '=================================='
echo -e "\nGet Fasta files\n"

cp $HOME/genomics/code/getall_fasta.py $TMPDIR

# arg1 sra_list location
# arg2 fasta path out
python3 getall_fasta.py $EPHEMERAL/sra_data/sra_runs.txt $EPHEMERAL/sra_data/files/

timer


#echo '=================================='
#echo -e "\nClear the environment\n"
#
#rm -f $HOME/ncbi/public/sra/*.lock # force for no error printing
#
#
#
#echo '=================================='
#echo -e "\nSelect the run code\n"
#
## array of data and select
#mapfile -t DATA < $DIR/sra_runs.txt 
##echo "${#DATA[@]}" # length of array
##FILE=${DATA[$PBS_ARRAY_INDEX]} # select the data by the job number
#
#for FILE in "${DATA[@]}"
#do
#
#
#echo '=================================='
#echo -e "\nFile: " $FILE
#echo -e '\n=================================='
#
#
#echo '=================================='
#echo -e "\nFetching\n"
#
## retrieve the SRA data in raw format
##prefetch $FILE --max-size 100G
#
#fastq-dump -X 5 -Z $FILE > $RES_DIR/$FILE'.fq'
#
#timer
#
#done
#echo '=================================='
#echo -e "\nMoving\n"
#
#mv $HOME/ncbi/public/sra/$FILE'.sra' $RES_DIR
#
#timer



#sam-dump $FILE | samtools view --threads 31 -bS - > $FILE.bam
#sam-dump ERR868003 | samtools view --threads 31 -bS - > $DIR/ERR868003.bam
#sam-dump ERR2179543 | samtools view --threads 31 -bS - > $DIR/ERR2179543.bam



