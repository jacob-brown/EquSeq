#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-15
# Last Modified: 2020-06-16
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

# plink conversion script - provided by treemix creators
wget https://bitbucket.org/nygcresearch/treemix/downloads/plink2treemix.py
