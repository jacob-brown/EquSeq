#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-23
# Last Modified: 2020-06-23



""" move wgs_data bam files to the final directory.
	some are in merged, and others sorted """

###########################################
################# Modules #################
###########################################
#module load anaconda3/personal # python
import csv
import subprocess

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

###########################################
######### Input(s) and Parameters #########
###########################################

no_merge = open_csv("data/no_merge.csv")
no_merge = open_csv("data/to_merge.csv")
sorted_files = subProWrap("ls /rds/general/user/jb1919/ephemeral/wgs_data/sorted/*")
merged_files = subProWrap("ls /rds/general/user/jb1919/ephemeral/wgs_data/merged/*")

###########################################
############### Wraggling #################
###########################################



###########################################
############### Analysis ##################
###########################################



###########################################
############### Plotting ##################
###########################################