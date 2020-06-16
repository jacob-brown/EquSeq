#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-15
# Last Modified: 2020-06-15
# Desc: 

cd dependancies
wget https://bitbucket.org/nygcresearch/treemix/downloads/treemix-1.13.tar.gz
tar -xvf treemix-1.13.tar.gz
cd treemix-1.13
./configure
make
make install
cd ../
rm treemix-1.13.tar.gz

# plink conversion
wget https://bitbucket.org/nygcresearch/treemix/downloads/plink2treemix.py

# tmp vcf from hpc
# on HPC 
cd $EPHEMERAL
head snp_calling/snps.chr3.raw.vcf -n50000 > sandbox/snps.vcf
# local
scp jb1919@login.cx1.hpc.ic.ac.uk:/rds/general/user/jb1919/ephemeral/sandbox/snps.vcf .

# plink - mac only
wget http://s3.amazonaws.com/plink1-assets/plink_mac_20200428.zip
mkdir plink
mv plink_mac_20200428.zip plink
cd plink
unzip plink_mac_20200428.zip
cd ../