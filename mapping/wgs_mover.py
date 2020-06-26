#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-23
# Last Modified: 2020-06-25


""" move wgs_data bam files to the final directory.
	some are in merged, and others sorted """

###########################################
################# Modules #################
###########################################


import csv
import subprocess
import os 
import sys

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

def subProWrap(command):
	""" wrapper for a subprocess command """
	p = subprocess.Popen([command], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	out, er = p.communicate()
	out_string = out.decode()
	files = out_string.replace("\t", ":").split("\n")
	files.remove("")
	return files

def listAbsFile(path):
	""" list full path of files using subProWrap() """
	tmp = subProWrap("ls {}*".format(path))
	#add_path = [path + i for i in tmp]
	return tmp #add_path

def stripPE(string):
	""" strip path and extension from 
		file string """
	noext = string.split(os.extsep, 1)[0]
	elem = os.path.split(noext)[1]
	return elem

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

no_merge = open_csv("/rds/general/user/jb1919/home/genomics/EquSeq/data/no_merge.csv")
no_runs = [i[1] for i in no_merge] # runs only

to_merge = open_csv("/rds/general/user/jb1919/home/genomics/EquSeq/data/to_merge.csv")
to_runs_nest = [i[1:] for i in to_merge] # runs only
to_runs = [item for sublist in to_runs_nest for item in sublist] # flattern 

# files in dirs
sorted_files = listAbsFile("ls /rds/general/user/jb1919/ephemeral/wgs_data/sorted/*")
merged_files = listAbsFile("ls /rds/general/user/jb1919/ephemeral/wgs_data/merged/*")

# lists from HPC
# sorted
#sorted_files = open_csv("sandbox/ls_all_sort.list")
#sorted_files = [i[0] for i in sorted_files] # flatten

# merged
#merged_files = open_csv("sandbox/ls_all_sort.list")
#merged_files = [i[0] for i in merged_files] # flatten

###########################################
############### Wraggling #################
###########################################

srt_bam = [i for i in sorted_files if ".bai" not in i]
srt_run_codes = [stripPE(i) for i in srt_bam]

# not in to be merged list
codes_from_sort = [i for i in srt_run_codes if i in no_runs]

# those in merge list
codes_from_merge = [i for i in srt_run_codes if i in to_runs]

# others - should be none
others = [i for i in srt_run_codes if i not in to_runs and i not in no_runs]


print("number of files missing from merge/sort search: " + str(len(others)))

if len(others) > 0:
	print("exiting.")
	sys.exit()

# get the path strings

# sorted
paths_from_sort = [i for i in sorted_files if stripPE(i) in codes_from_sort]

if len(paths_from_sort)/2 == len(codes_from_sort): 
	print("path list correct length")
else:
	print("path and code list are not the same length")
	sys.exit()

# merged - should be fine as is
#merged_files

# combine 
final_paths = paths_from_sort + merged_files

# save list
#files = listAbsFile("sandbox/moving/fol1/")
saveTxt("/rds/general/user/jb1919/ephemeral/wgs_data/move.list", final_paths)






