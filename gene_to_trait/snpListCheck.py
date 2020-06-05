#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-02
# Last Modified: 2020-06-02



""" Check which gene variants are present in the snp lists.
	snp lists generated after samtools mpileup """

###########################################
################# Modules #################
###########################################

import csv
import subprocess as sp

###########################################
############## Function(s) ################
###########################################

def open_csv(file):
	
	""" open a csv into a list format """

	tmp = [] # initialise the list
	with open(file, 'r') as f:
		reader = csv.reader(f)
		for row in reader:
			tmp.append(row) # add row to list

	return tmp


###########################################
######### Input(s) and Parameters #########
###########################################

# varID, chr, phen, breed, position, 
	# position_stop, type, reference_nucleotide, new_nucleotide 
snp_info = open_csv("data/gene_variants/snp.list.info.csv")

command = "cat data/processed_sequences/snp.list.dump/*"
p = sp.Popen([command], stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
out, er = p.communicate()
out_decode = out.decode()
pos_list = out_decode.split("\n")



###########################################
############### Wraggling #################
###########################################

to_match = ["chr" + i[1] + "_" + i[4] for i in snp_info]

# none exist...
[i for i in pos_list if i in to_match]


###########################################
############### Analysis ##################
###########################################



###########################################
############### Plotting ##################
###########################################