#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-03-03
# Last Modified: 2020-03-03



""" Match pair ended files and run alignment using BWA """

###########################################
################# Modules #################
###########################################


import os 
import re
import sys

###########################################
######### Input(s) and Parameters #########
###########################################

# path to trimmed reads
path_in = '/rds/general/user/jb1919/ephemeral/trimmed'

# path to output files
path_out = '/rds/general/user/jb1919/ephemeral/aligned'

# path to indexed reference genome
ref_genome = '/rds/general/user/jb1919/ephemeral/ref_genome_indexed/EquCab2.fna'

###########################################
############ File Preparation #############
###########################################


# list of files
files = os.listdir(path_in)  
files.sort()

# match everything after AAAAA and before .fq
pattern_read = r'.*(?=_\d.fq)'# all the read info



### strip the file names ###
store = []

for i in files:
	
	# extract all info - read and pair
	match_read = re.search(pattern_read, i)
	read = match_read.group()

	# update store
	store.append([i, read])

### match the pairs ###
pairs = []

for c in store:

	# match the pairs
	p = [i[0] for i in store if i[1] == c[1]]
	
	# file out name - SAM 
	p_out = '{}.sam'.format(c[1]) 

	# append
	p.append(p_out) # to list 
	
	# only append when not in th elist already
		# matches will return doubles
	if p not in pairs:
		pairs.append(p) # to storage


###########################################
################ Command ##################
###########################################

# catch the job number
job_num = int(os.environ['PBS_ARRAY_INDEX'])

selected = pairs[job_num]

### full paths ### 
read1 = '{}/{}'.format(path_in, selected[0])
read2 = '{}/{}'.format(path_in, selected[1])
out_file = '{}/{}'.format(path_out, selected[2])

# construc the command
command = 'bwa mem {} {} {} > {}'.format(ref_genome, read1, read2, out_file) 
# bwa mem ref.fa read1.fq read2.fq > aln-pe.sam


### run alignment ###
os.system(command)

