#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-09
# Last Modified: 2020-06-09
# Desc: 


# generate merged log file

loc=$1


(for log in `ls $loc/*.log`;
do 
	K=($(ggrep -Po 'ALL.MIX.K\K[^ ]+' $log));
	VAL=($(ggrep -Po 'like=\K[^ ]+' $log));
	echo $K"\t"$VAL
done) > $loc/logfile.merge



















# guide
# https://github.com/alexkrohn/AmargosaVoleTutorials/blob/master/ngsAdmix_tutorial.md
#cd sandbox
#mkdir logfiles
#cp ../results/ancestry/eu_more_wg_5kb_05maf/*.log logfiles/
#
#cd logfiles
## use ggrep rather than grep - mac issues?
#(for log in `ls *.log`; do ggrep -Po 'like=\K[^ ]+' $log; done) > ../logfile
#cd ../
#
## R
#R
#logs <- as.data.frame(read.table("logfile"))
#
## number of K values used should reflect the reps
#	# we used 5 so...
#logs$K <- c(rep("2", 1), rep("3", 1), rep("4", 1), rep("5", 1), rep("6", 1))
#write.table(logs[, c(2, 1)], "logfile_formatted", row.names = F, col.names = F, quote = F)
#
#
## bash
## install clump
#wget http://clumpak.tau.ac.il/download/CLUMPAK.zip
#unzip CLUMPAK.zip
##cd CLUMPAK
#
##cd ../
#BESTK=CLUMPAK/26_03_2015_CLUMPAK/CLUMPAK/BestKByEvanno.pl
#
#cd CLUMPAK/26_03_2015_CLUMPAK/CLUMPAK/
#perl BestKByEvanno.pl --i 1 --d out --f ../../../logfile_formatted -- inputtype lnprobbyk#