#!/bin/bash
#PBS -lwalltime=05:00:00
#PBS -lselect=1:ncpus=24:ompthreads=24:mem=240gb


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
module load fastx/0.0.14
source activate myenv # activate conda environment


DBNAME=$EPHEMERAL/kraken/db_kraken_horse


echo '=================================='
echo -e "\nBuild taxonomy database\n"

# pull taxonomy structure 
kraken2-build --download-taxonomy --db $DBNAME

echo '=================================='
echo -e "\nAdd libraries\n"
# add reference libraries

echo -e "\narchaea\n"
kraken2-build --download-library archaea --db $DBNAME
echo -e "\nbacteria\n"
kraken2-build --download-library bacteria --db $DBNAME
echo -e "\nviral\n"
kraken2-build --download-library viral --db $DBNAME
echo -e "\nfungi\n"
kraken2-build --download-library fungi --db $DBNAME
echo -e "\nplant\n"
kraken2-build --download-library plant --db $DBNAME
#echo -e "\nhuman\n"
#kraken2-build --download-library human --db $DBNAME

# update taxanomic info in fasta file
# taxanomic IDs from:
	#https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi
# Update header 
#var=">EquCab2|kraken:taxid|9796  Horse reference genome"
# line 1 and change
#sed "1c$var" $EPHEMERAL/kraken/EquCab2.fna > horse.fa

#echo -e "\nhorse\n"
kraken2-build --add-to-library $EPHEMERAL/kraken/horse.fa --db $DBNAME

echo '=================================='
echo -e "\nFinal Build\n"

# build the database
kraken2-build --build --threads 24 --db $DBNAME

echo '=================================='
echo -e "\nclose conda\n"

conda deactivate
