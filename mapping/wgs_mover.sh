#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-24
# Last Modified: 2020-06-25
# Desc: move files in move.list to final directory

#module load anaconda3/personal 

FILES=($(cat $EPHEMERAL/wgs_data/move.list))
OUT_DIR=$EPHEMERAL/wgs_data/final

for file in "${FILES[*]}"
do
	#echo "moving " $file " to " $OUT_DIR "\n"
	mv $file $OUT_DIR
done

