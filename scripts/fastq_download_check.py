#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-12
# Last Modified: 2020-06-12



"""  """

###########################################
################# Modules #################
###########################################

import csv

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

f = open("sandbox/fasta_logs/log_fail/issue.run.codes", "r")
data=f.readlines()
issue_codes = [i.replace("\n", "") for i in data]

info_all = open_csv("data/cleaned_data/info_all.csv")

###########################################
############### Wraggling #################
###########################################

issue_info = [i for i in info_all if i[0] in issue_codes]

write_csv(issue_info, "sandbox/fasta_logs/issue_info.csv")

print("saving: data/cleaned_data/rerun.codes.txt")
print("first delete the files on the HPC, then rerun these codes.")

saveTxt("data/cleaned_data/rerun.codes.txt", issue_codes)













