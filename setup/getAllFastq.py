#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-04-16
# Last Modified: 2020-06-11



""" Which fasta files do we still need? """

###########################################
################# Modules #################
###########################################

import csv
import os
import numpy as np

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

def write_csv(list_file, path):
	
	""" Write list to csv """

	with open(path, 'w') as f:
		writer = csv.writer(f, delimiter=',')
		for i in list_file:
			writer.writerow(i)

###########################################
######### Input(s) and Parameters #########
###########################################

# if sra_runs is present filter
if os.path.isfile('data/cleaned_data/sra_runs.txt'):
	code_file = 'data/cleaned_data/sra_runs.txt'

	with open(code_file, 'r') as f:
		codes = f.readlines() # read file
else:
	codes = []


info_all_path = 'data/cleaned_data/info_all.csv'

info_all = open_csv(info_all_path)

###########################################
############### Wraggling #################
###########################################


# only modern samples
i_era = info_all[0].index('era')
modern = [i for i in info_all if i[i_era] == 'modern']

# samples we don't have
# codes in sra file (those we have)
codes = [i.strip('\n') for i in codes] # filter and strip new line

i_run = info_all[0].index('Run')
i_org = info_all[0].index('Organism')
i_sub = info_all[0].index('sub_group')

#modern_unused = []
#for elm in modern:
#	if elm[i_run] not in codes:
#		if elm[i_sub] == '' and elm[i_org] == 'Equus caballus':
#			# not a generic horse 
#			pass
#		else:
#			modern_unused.append(elm)


# not already in code list
modern_unused = [i for i in modern if i[i_run] not in codes] 

# not a generic horse - remove ungrouped
modern_unused = [i for i in modern_unused if i[i_sub]!= '']

modern_unused.insert(0, info_all[0])

# codes only
codes_to_use = [i[0] for i in modern_unused[1:]]


# all data still to use
write_csv(modern_unused, "data/cleaned_data/breeds_to_use.csv")

# codes only
np.savetxt('data/cleaned_data/sra_runs_to_do.txt', codes_to_use, fmt='%s')


