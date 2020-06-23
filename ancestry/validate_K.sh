#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-09
# Last Modified: 2020-06-23
# Desc: used within validate_K.R


# generate merged log file
	# ggrep usage on MacOS - grep for gnu

# log files dir
loc=$1 

(for log in `ls $loc/*.log`;
do 
	K=($(ggrep -Po 'ALL.MIX.K\K[^ ]+' $log));
	VAL=($(ggrep -Po 'like=\K[^ ]+' $log));
	echo $K"\t"$VAL
done) > $loc/logfile.merge
