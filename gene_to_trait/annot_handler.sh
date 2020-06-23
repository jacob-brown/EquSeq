#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-18
# Last Modified: 2020-06-18
# Desc: combine and anylise vep output


for i in {0..9}
do
	cat annot.$i.raw.txt | sed -n '/#Uploaded_variation/,$p' >> combined.txt
done
