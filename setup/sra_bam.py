#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-04-09
# Last Modified: 2020-04-09


###### REMOVE this file #######

""" check if bam files are present  """

###########################################
################# Modules #################
###########################################

import getall_fasta 


###########################################
############## Function(s) ################
###########################################



###########################################
######### Input(s) and Parameters #########
###########################################

files = '../data/cleaned_data/sra_runs.txt'

with open(files) as f:
    runs = f.read().splitlines()

###########################################
############### Wraggling #################
###########################################

store = []
for i in runs:
	tmp = getall_fasta.getTableInfo(i, 'submitted_format', 'submitted_ftp')
	tmp.insert(0, i)
	store.append(tmp)


store[10]

# bams
poss = [i for i in store if len(i) == 3]
bam = [i for i in poss if 'BAM' in i[2][0].upper()]

# no bams
no_bam = [i for i in poss if 'BAM' not in i[2][0].upper()]
nob = [i for i in store if len(i) != 3]
nob.append(no_bam)


###########################################
############### Analysis ##################
###########################################



###########################################
############### Plotting ##################
###########################################