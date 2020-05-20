#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-05-20
# Last Modified: 2020-05-20
# Desc: install taxonkit, one time run


cd data/oral_diversity/
brew install brewsci/bio/taxonkit # for mac
taxonkit genautocomplete

wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

tar xvzf taxdump.tar.gz

mkdir $HOME/.taxonkit
cp names.dmp nodes.dmp delnodes.dmp merged.dmp $HOME/.taxonkit