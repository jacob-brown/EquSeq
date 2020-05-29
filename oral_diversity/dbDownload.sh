#!/bin/bash
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=1:mem=1gb

# download data into kraken2 database

source $HOME/genomics/EquSeq/scripts/unix_functions.sh

echo '=================================='
echo -e "\nLoad modules\n"

module load anaconda3/personal
source activate myenv 


DIR=$EPHEMERAL/oral_diversity/
DBNAME=$DIR/db_kraken_horse

echo '=================================='
echo -e "\nBuild taxonomy database\n"

# pull taxonomy structure 
kraken2-build --download-taxonomy --db $DBNAME

timer

echo '=================================='
echo -e "\nAdd libraries\n"

echo -e "\narchaea\n"
kraken2-build --download-library archaea --db $DBNAME

timer

echo -e "\nbacteria\n"
kraken2-build --download-library bacteria --db $DBNAME

timer

echo -e "\nviral\n"
kraken2-build --download-library viral --db $DBNAME

timer

echo -e "\nfungi\n"
kraken2-build --download-library fungi --db $DBNAME

timer

echo -e "\nplant\n"
kraken2-build --download-library plant --db $DBNAME

timer

echo -e "\nhorse\n"

# add to database
kraken2-build --add-to-library $DIR/horse.fa --db $DBNAME

timer

echo -e "\nhuman\n"
kraken2-build --add-to-library $DIR/homo_sap.fa --db $DBNAME

timer

conda deactivate





