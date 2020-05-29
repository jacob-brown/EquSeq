#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-05-25
# Last Modified: 2020-05-27
# Desc: download horse and human genomes, preparing them for kraken db


DIR=$EPHEMERAL/oral_diversity/

echo '=================================='
echo -e "\nhuman ref genome\n"

rsync --copy-links --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz $DIR

gunzip $DIR/*genomic.fna.gz 

echo "update header"
var=">GRCh38.p13|kraken:taxid|9606  Human reference genome" \
	&& sed "1c$var" $DIR/GCF_000001405.39_GRCh38.p13_genomic.fna > homo_sap.fa


echo '=================================='
echo -e "\nhorse ref genome\n"

cp $EPHEMERAL/ref_genome/EquCab3.fna $DIR

echo "update header"
var=">EquCab3|kraken:taxid|9796  Horse reference genome" \
	&& sed "1c$var" $DIR/EquCab3.fna > horse.fa

echo "tidy up"
rm -f $DIR/EquCab3.fna $DIR/GCF_000001405.39_GRCh38.p13_genomic.fna