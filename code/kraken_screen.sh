#!/bin/bash
#PBS -lwalltime=00:30:00
#PBS -lselect=1:ncpus=1:mem=9gb

# 0-5 doing 
#15-23 done - under wrong files

# wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken2_v1_8GB_201904_UPDATE.tgz

#wget http://github.com/DerrickWood/kraken2/archive/v2.0.8-beta.tar.gz

# gunzip minikraken2_v1_8GB_201904_UPDATE.tgz
# gunzip v2.0.8-beta.tar.gz

#tar -xvf minikraken2_v1_8GB_201904_UPDATE.tar
#tar -xvf v2.0.8-beta.tar

DATA=($EPHEMERAL/aligned/*) # array of all data
FILE=${DATA[$PBS_ARRAY_INDEX]} # select the data by the job number
NOEXT=$(echo $FILE | cut -f 1 -d '.') # remove extension 
BASE=$(basename "$NOEXT") # remove path

FILE_1=$EPHEMERAL/reads/$BASE'_1.fq'
FILE_2=$EPHEMERAL/reads/$BASE'_2.fq'

echo $FILE_1
echo $FILE_2

#----- load modules ----#
echo '=================================='
echo -e "\nLoad modules\n"
module load anaconda3/personal
source activate myenv # activate conda environment



DIR=$EPHEMERAL/kraken/

echo '=================================='
echo -e "\nRun Kraken\n"


#kraken2 --db $DIR/db_kraken_horse \
#		--threads 24 \
#		--report $DIR/MYSAMPLE.full.report \
#		--paired $DIR/V300044309_L2_B5GHORlfyRAAAAAAA-517_1.fq \
#		$DIR/V300044309_L2_B5GHORlfyRAAAAAAA-517_2.fq \
#		> $DIR/MYSAMPLE.full.kraken

		#$DIR/V300044309_L2_B5GHORlfyRAAAAAAA-517_1.fq \
		#$DIR/V300044309_L2_B5GHORlfyRAAAAAAA-517_2.fq \

kraken2 --db $DIR/minikraken2_v2_8GB \
		--report $DIR/kreport_out/$BASE.mini.report \
		--quick \
		--fastq-input \
		--paired $FILE_1 $FILE_2 \
		> $DIR/kraken_out/$BASE.mini.kraken


# close environment
conda deactivate



# pavian  for quick view
#https://ccb.jhu.edu/software/pavian/index.shtml
# R
#pavian::runApp()