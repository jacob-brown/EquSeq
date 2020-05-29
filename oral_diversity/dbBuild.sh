#!/bin/bash
#PBS -lwalltime=05:00:00
#PBS -lselect=1:ncpus=20:ompthreads=20:mem=240gb

echo '=================================='
echo -e "\nLoad modules\n"
module load anaconda3/personal
source activate myenv # activate conda environment

DIR=$EPHEMERAL/oral_diversity/
DBNAME=$DIR/db_kraken_horse


echo '=================================='
echo -e "\nFinal Build\n"

# build the database
kraken2-build --build --threads 19 --db $DBNAME

echo '=================================='
echo -e "\nclose conda\n"

conda deactivate
