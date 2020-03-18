#!/bin/bash
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=1:mem=1gb

# download data into kraken2 database


#----- load modules ----#
echo '=================================='
echo -e "\nLoad modules\n"
module load anaconda3/personal
module load fastx/0.0.14
source activate myenv # activate conda environment


DBNAME=$EPHEMERAL/kraken/db_kraken_horse


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

echo -e "\nhorse\n"
# update taxanomic info in fasta file
# taxanomic IDs from:
	#https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi
# Update header 
#var=">EquCab2|kraken:taxid|9796  Horse reference genome"
#sed "1c$var" $EPHEMERAL/kraken/EquCab2.fna > horse.fa # line 1 and change
# add to database
kraken2-build --add-to-library $EPHEMERAL/kraken/horse.fa --db $DBNAME

echo -e "\nhuman\n"
kraken2-build --add-to-library $EPHEMERAL/kraken/homo_sap.fa --db $DBNAME

echo '=================================='
echo -e "\nclose conda\n"

conda deactivate





