#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-03-02
# Last Modified: 2020-03-02
# Desc: 


# EquCab2 
# ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Equus_caballus/all_assembly_versions/suppressed/

# filename: GCF_000002305.2_EquCab2.0_genomic.fna.gz

# EquCab3 location:
# ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Equus_caballus/latest_assembly_versions/

#-------------- EquCab2 --------------#

# specific file 
rsync --copy-links --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Equus_caballus/all_assembly_versions/suppressed/GCF_000002305.2_EquCab2.0/GCF_000002305.2_EquCab2.0_genomic.fna.gz .

# all files
#rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Equus_caballus/all_assembly_versions/suppressed/GCF_000002305.2_EquCab2.0 .

#-------------- EquCab3 --------------#
# specific file 
rsync --copy-links --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Equus_caballus/latest_assembly_versions/GCF_002863925.1_EquCab3.0/GCF_002863925.1_EquCab3.0_genomic.fna.gz .

#-------------- unzip --------------#

gunzip GCF_000002305.2_EquCab2.0_genomic.fna.gz

#-------------- rename --------------#
mv GCF_000002305.2_EquCab2.0_genomic.fna EquCab2.fna
#mv GCF_002863925.1_EquCab3.0_genomic.fna EquCab3.fna