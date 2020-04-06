#!/bin/bash
#PBS -lwalltime=02:00:00
#PBS -lselect=1:ncpus=10:ompthreads=10:mem=120gb

# ~30 min per pair
# 0-23 jobs

# use aligned as only _1 and _2 extensions are absent
DATA=($EPHEMERAL/aligned/*) # array of all data
FILE=${DATA[$PBS_ARRAY_INDEX]} # select the data by the job number
NOEXT=$(echo $FILE | cut -f 1 -d '.') # remove extension 
BASE=$(basename "$NOEXT") # remove path

# files 
FILE_1=$EPHEMERAL/reads/$BASE'_1.fq'
FILE_2=$EPHEMERAL/reads/$BASE'_2.fq'

# database
DBNAME=$EPHEMERAL/kraken/db_kraken_horse/

# directory
DIR=$EPHEMERAL/kraken/

#----- load modules ----#
echo '=================================='
echo -e "\nLoad modules\n"
module load anaconda3/personal
source activate myenv # activate conda environment

DIR=$EPHEMERAL/kraken/

echo '=================================='
echo -e "\nRun Kraken\n"


kraken2 --db $DBNAME \
		--report $DIR/kreport_out/$BASE.kreport \
		--fastq-input \
		--threads 10 \
		--paired $FILE_1 $FILE_2 \
		> $DIR/kraken_out/$BASE.kraken

echo '=================================='
echo -e "\nclose conda\n"

# close environment
conda deactivate



# pavian  for quick view
#https://ccb.jhu.edu/software/pavian/index.shtml
# R
#pavian::runApp()