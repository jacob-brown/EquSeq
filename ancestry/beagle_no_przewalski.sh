#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-07-06
# Last Modified: 2020-07-13
# Desc: Remove Przewalski (and hybrids) from beagle file and run a PCA



# problem ids
# Ind84 
# Ind85 
# Ind86 
# Ind87 
# Ind88
# Ind89 
# Ind114
# Ind115
# Ind118
# Ind119


cp ALL.merged.beagle.gz ../sandbox
cd ../sandbox

zcat ALL.merged.beagle.gz | head -n1 > tmp.header


# 256-273 are the problem columns
cat tmp.header | tr '\t' '\n' | grep -n Ind8

# 346-351, 358-363
cat tmp.header | tr '\t' '\n' | grep -n Ind11
cat tmp.header | tr '\t' '\n' | wc -l

# look and reslove Ind84-89
cat tmp.header | cut -f1-255,274-519 | grep Ind8

# look and reslove Ind114-Ind115 
cat tmp.header | cut -f1-345,352-519 | grep Ind11

# look and reslove Ind118-Ind119
cat tmp.header | cut -f1-357,364-519 | grep Ind11


# look and reslove all combined 
cat tmp.header | cut  -f1-255,274-345,352-357,364-519 | grep Ind11
cat tmp.header | cut  -f1-255,274-345,352-357,364-519 | grep Ind8

# if all good continue


zcat ALL.merged.beagle.gz | cut -f1-255,274-345,352-357,364-519 > NO.PRZ.merged.beagle 
gzip NO.PRZ.merged.beagle
mv NO.PRZ.merged.beagle.gz ../ancestry/


