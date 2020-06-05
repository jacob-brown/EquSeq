#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-03
# Last Modified: 2020-06-05
# Desc: 

### setup in conda ###
module load anaconda3/personal

conda install -c bioconda ensembl-vep

mkdir $HOME/.vep
cd $HOME/.vep
wget ftp://ftp.ensembl.org/pub/release-100/variation/indexed_vep_cache/equus_caballus_vep_100_EquCab3.0.tar.gz
tar xzf equus_caballus_vep_100_EquCab3.0.tar.gz

mkdir ~/.vep/Plugins/