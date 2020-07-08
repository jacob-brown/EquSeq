#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-07-06
# Last Modified: 2020-07-06


""" prepare snp list generated from initial bam 
run for use in generating vcf files """



###########################################
################# Modules #################
###########################################

from natsort import natsorted
import numpy as np

###########################################
############## Function(s) ################
###########################################

### save a text file without a new line at the end
def saveTxt(dirfile, listToSave, sep='\n'):
	with open(dirfile, 'w') as f:
		for num, val in enumerate(listToSave):
			if num == len(listToSave) - 1:
				f.write(val)
			else:
				f.write(val + sep)


###########################################
######### Input(s) and Parameters #########
###########################################

fopen = open("data/ancestry/snp.all.list", "r")
fdata = fopen.readlines()
snplist =[i.replace("\n", "") for i in fdata]
n = 50

###########################################
############### Wraggling #################
###########################################

# sort
snplist_sorted = natsorted(snplist, key=lambda y: y)
snplist_sorted = [i.replace("_", ":") for i in snplist_sorted]

# break down into n files
split_snps = np.array_split(np.array(snplist_sorted), 50)

# save 
for num, val in enumerate(split_snps):
	fileOut = "data/ancestry/snp.vcf.list/snp.{}.vcf.list".format(num)
	saveTxt(fileOut, val)







