#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-19
# Last Modified: 2020-06-19
# Desc: 



DIR=sandbox/beagle_traits/
beagle=data/processed_sequences/beagle/ALL.merged.beagle.gz
SNP=data/gene_variants/trait.snps/trait.snp.mend.list

### all inds ###
head $SNP

zcat < $beagle | grep chr1_95257460
zcat < $beagle | grep -Fw -f $SNP > $DIR/tmp.beagle


zcat < $beagle | cut -f1 | grep -Fw -f $SNP > $DIR/tmp.match

zcat < $beagle | cut -f1 > $DIR/beagle.snps.txt

### benson only ###
zcat < $beagle | head -n1 
zcat < $beagle | cut -f1-3,136-138 | gzip > $DIR/benson.beagle.gz

zcat < $DIR/benson.beagle.gz | grep chr13_20501489

snparray=($(cat $SNP))
echo "${snparray[0]}"

for i in "${snparray[*]}"
do
	zcat < $DIR/benson.beagle.gz | grep $i 
done


python3 ancestry/snpHandler.py -c subBeagle -i $DIR/benson.beagle.gz -o $DIR/trait.snps -l $SNP

head sandbox/beagle_traits/trait.snps.beagle

"zcat {} | head -1 > {}.TMP.HEAD && zcat {} | grep -Fw -f {} - | cat {}.TMP.HEAD - > {} ".format(beagleIn, beagleIn, beagleIn, posList, beagleIn, outFile)

zcat < $DIR/benson.beagle.gz | grep -Fw -f $SNP > $DIR/tmp.beagle


head $DIR/tmp.beagle
zcat < $DIR/benson.beagle.gz | cut -f4 | sort | uniq