#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-08-06
# Last Modified: 2020-08-06
# Desc: 


samtools depth final.bam -r chr2 > ../../../sandbox/depthout2.csv
samtools coverage final.bam > ../../../sandbox/coverage.csv