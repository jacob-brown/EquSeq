#!/bin/bash
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=7:mem=10gb

# wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken2_v1_8GB_201904_UPDATE.tgz

#wget http://github.com/DerrickWood/kraken2/archive/v2.0.8-beta.tar.gz

# gunzip minikraken2_v1_8GB_201904_UPDATE.tgz
# gunzip v2.0.8-beta.tar.gz

#tar -xvf minikraken2_v1_8GB_201904_UPDATE.tar
#tar -xvf v2.0.8-beta.tar


#----- load modules ----#
echo '=================================='
echo -e "\nLoad modules\n"
module load anaconda3/personal
source activate myenv # activate conda environment



DIR=$EPHEMERAL/kraken/

echo '=================================='
echo -e "\nRun Kraken\n"

kraken2 --db $DIR/minikraken2_v1_8GB \
		--threads 6 \
		--report $DIR/MYSAMPLE.kreport2 \
		--paired $DIR/V300044309_L2_B5GHORlfyRAAAAAAA-517_1.fq \
		$DIR/V300044309_L2_B5GHORlfyRAAAAAAA-517_2.fq \
		> $DIR/MYSAMPLE.kraken2 


# close environment
conda deactivate


#jb1919@login.cx1.hpc.ic.ac.uk:/rds/general/user/jb1919/ephemeral/kraken


# KRONA

# pavian  for quick view
#https://ccb.jhu.edu/software/pavian/index.shtml
# R
#pavian::runApp(port=5000)