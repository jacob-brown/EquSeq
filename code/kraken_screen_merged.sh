#!/bin/bash
#PBS -lwalltime=70:00:00
#PBS -lselect=1:ncpus=10:ompthreads=10:mem=120gb

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


# database
DBNAME=$EPHEMERAL/kraken/db_kraken_horse/

# directory
DIR=$EPHEMERAL/kraken/

#----- load modules ----#
echo '=================================='
echo -e "\nLoad modules\n"
module load anaconda3/personal
module load java/jdk-8u144 # picard associated
module load picard/2.6.0 # cleaning

PICARD=$PICARD_HOME/picard.jar
source activate myenv # activate conda environment


FILE=$DIR/file.fa

timer

echo '=================================='
echo -e "\nconvert ot fasta\n"

java -Xmx120g -jar $PICARD SamToFastq \
		INPUT=$EPHEMERAL/mapping/merged/new.rg.bam \
		FASTQ=$DIR/  \
		TMP_DIR=$TMPDIR # resolves memory issues

timer

echo '=================================='
echo -e "\nRun Kraken\n"


kraken2 --db $DBNAME $FILE \
		--report $FILE.kreport \
		--fastq-input \
		--threads 10 > \
		$FILE.kraken

timer

echo '=================================='
echo -e "\nclose conda\n"

# close environment
conda deactivate



# pavian  for quick view
#https://ccb.jhu.edu/software/pavian/index.shtml
# R
#pavian::runApp()