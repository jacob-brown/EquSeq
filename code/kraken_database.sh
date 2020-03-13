#!/bin/bash
#PBS -l walltime=03:00:00
#PBS -l select=1:ncpus=8:mem=10gb



#install and create in a segregated conda environment. i.e.
#$ module load anaconda3/personal
#$ conda create -n myenv kraken2 
	# create myenv an install kraken2
#$conda install --name myenv libiconv
	# add package libiconv
#$ source activate myenv



#----- load modules ----#
echo '=================================='
echo -e "\nLoad modules\n"
module load anaconda3/personal
source activate myenv # activate conda environment


DBNAME=$EPHEMERAL/kraken/db_kraken_horse

# create a kraken database

echo '=================================='
echo -e "\nBuild taxonomy database\n"

# pull taxonomy structure 
kraken2-build --download-taxonomy --db $DBNAME

echo '=================================='
echo -e "\nAdd libraries\n"
# add reference libraries

echo -e "\narchaea\n"
#kraken2-build --download-library archaea --db $DBNAME
echo -e "\nbacteria\n"
kraken2-build --download-library bacteria --db $DBNAME
echo -e "\nviral\n"
#kraken2-build --download-library viral --db $DBNAME
echo -e "\nfungi\n"
#kraken2-build --download-library fungi --db $DBNAME
echo -e "\nplant\n"
#kraken2-build --download-library plant --db $DBNAME
echo -e "\nhuman\n"
#kraken2-build --download-library human --db $DBNAME
echo -e "\nhorse\n"
kraken2-build --download-library horse --db $DBNAME

echo '=================================='
echo -e "\nFinal Build\n"

# build the database
kraken2-build --build --threads 7 --db $DBNAME

echo '=================================='
echo -e "\nclose conda\n"

conda deactivate
