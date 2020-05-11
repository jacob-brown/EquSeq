#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-05-07
# Last Modified: 2020-05-07
# Desc: 


DIR=$1

# unzip all 
gunzip $DIR/*.beagle.gz
echo "unzipped"
FILES=($( ls -v $DIR/*.beagle ))

echo "combining beagle files"
head -1 ${FILES[0]} > $DIR/tmp.head
awk FNR-1 $DIR/*.beagle | cat $DIR/tmp.head - \
	> $DIR/ALL.merged.beagle && rm $DIR/tmp.head

# zip again
gzip $DIR/*.beagle
echo "zipped and complete"