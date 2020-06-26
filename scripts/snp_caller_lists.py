#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-23
# Last Modified: 2020-06-25



""" generate a multiple snp lists for samtools to 
	use for mpileup. Quicker than per chromosome runs. """


# first run: snp_caller_lists.sh

###########################################
################# Modules #################
###########################################

import numpy as np
import math

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

fopen = open("data/EquCab3.stats.txt")
file_str = fopen.readlines()
data = [i.replace("\n", "").split("\t") for i in file_str]
# only chromosomes and their length, no X chromo
data_chrom = [i for i in data 
				if i[1] != 'X'
					and i[2] == 'Chromosome' 
					and i[4] == 'total-length' 
					and i[3] == 'all']

des_len = 10000000 # desired number of snps per file

###########################################
############### Wraggling #################
###########################################
 
lengths = [int(i[5]) for i in data_chrom]

# split into multiple files (.bed)
	# chr1 0 100 - start is inclusive, end is non inclusive

store = []
for elem in data_chrom:

	chrom = "chr" + elem[1]
	length = int(elem[5])
	fileN = math.ceil(length/des_len) # number of files needed
	start_stop = np.linspace(0, length, dtype = int, num=fileN)
	start_stop = start_stop.astype('str')
	start = start_stop[:-1]
	stop = start_stop[1:]

	for i in range(0, len(start)):
		string = chrom + "\t" + start[i] + "\t" + stop[i]
		store.append(string)
		path = "data/snp_calling_list/snpcall.{}.{}.list.bed".format(chrom, i)
		wfile = open(path, 'w')
		wfile.writelines(string)
		

# rm data/snp_calling_lists/*


###########################################
############### Analysis ##################
###########################################



###########################################
############### Plotting ##################
###########################################