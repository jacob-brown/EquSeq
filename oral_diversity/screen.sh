#!/bin/bash
#PBS -lwalltime=02:00:00
#PBS -lselect=1:ncpus=10:ompthreads=10:mem=120gb

# ~30 min per pair
# 0-23 jobs

#qsub -J 0-23 screen.sh
#qsub -J 0-2 screen.sh

# use aligned as only _1 and _2 extensions are absent
# use combined pairs as a ref for the base file names
DATA=($EPHEMERAL/novel_data/converted/*) 
FILE=${DATA[$PBS_ARRAY_INDEX]} # select the data by the job number
NOEXT=$(echo $FILE | cut -f 1 -d '.') # remove extension 
BASE=$(basename "$NOEXT") # remove path

# files 
FILE_1=$EPHEMERAL/novel_data/raw_files/$BASE'_1.fq'
FILE_2=$EPHEMERAL/novel_data/raw_files/$BASE'_2.fq'


echo "file 1: " $FILE_1
echo "file 2: " $FILE_2

# database
DBNAME=$EPHEMERAL/oral_diversity/db_kraken_horse/

# directory
DIR=$EPHEMERAL/oral_diversity/

#----- load modules ----#
echo '=================================='
echo -e "\nLoad modules\n"
module load anaconda3/personal
source activate myenv # activate conda environment


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

### braken ###
#cd $EPHEMERAL/sandbox/
##conda install -c bioconda bracken
#git clone https://github.com/jenniferlu717/Bracken.git
#cd Bracken
#bash install_bracken.sh
#cd ..
#
## testing
#scp V300044309_L2_B5GHORlfyRAAAAAAA-517.kreport jb1919@login.cx1.hpc.ic.ac.uk:/rds/general/user/jb1919/ephemeral/sandbox
## S for species level
#Bracken/bracken -d ../oral_diversity/db_kraken_horse -i V300044309_L2_B5GHORlfyRAAAAAAA-517.kreport -o bracken.species.txt -l S
#






