#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-23
# Last Modified: 2020-06-23
# Desc: 

FILE_IN=data/EquCab3.0_assembly_stats.txt

wget -O $FILE_IN https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Equus_caballus/all_assembly_versions/GCF_002863925.1_EquCab3.0/GCF_002863925.1_EquCab3.0_assembly_stats.txt


sed '1,/unit-name/d' $FILE_IN > data/EquCab3.stats.txt
