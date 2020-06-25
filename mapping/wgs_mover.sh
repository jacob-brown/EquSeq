#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-24
# Last Modified: 2020-06-24
# Desc: 


FILES=($(cat sandbox/moving/move.list))
OUT_DIR=sandbox/moving/fol2/

for file in "${FILES[*]}"
do
	echo "moving " $file " to " $OUT_DIR "\n"
	#mv $file $OUT_DIR
done
