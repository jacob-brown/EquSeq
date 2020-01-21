#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-01-14
# Last Modified: 2020-01-20



""" Gather SRA infotable for each run """

###########################################
################# Modules #################
###########################################

import subprocess
import numpy as np
import csv
import re

###########################################
############## Function(s) ################
###########################################


def split_comma(string):

	""" commas within text comments caused issues with the split
		this function splits on index and removes commas """

	p = re.compile(r',(?=([^"]*"[^"]*")*[^"]*$)')

	# determine index string shoule be broken down in
	index_break = [m.start() for m in p.finditer(string)] 
	index_break.insert(0,0) # first index should be 0, note the -1 in loop below

	store_elem = [] # store the results

	for i in range(1,len(index_break)):
		p1 = index_break[i-1] # position 1 
		p2 = index_break[i] # position 2
		strg = string[p1:p2].replace(",","") # remove all commas
		store_elem.append(strg) # store the data

	return store_elem


def get_runinfo(code, header_only=False):
	
	""" gather infotables on run sequences using entrez query. 
		previous issues with the header resulted in seperation of the two into header and body """

	command = "esearch -db sra -query {} | efetch -format runinfo".format(code) # entrez query to send to terminal
	p = subprocess.Popen(command, shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) # run command and store the results
	stdout, stderr = p.communicate() 

	# blank returns
	if stdout.decode() != '':

		tmp = stdout.decode().split('\n') # each infotable is at a linebreak
		tmp = [i for i in tmp if len(i) != 0] # remove blank entries generated from linebreaks

		header = tmp[0]
		body = tmp[1:]

		### break down into a list of lists ###
		store = [] 

		for elem in body: 

			if elem != header: # prevent headers from getting generated in data body

			    body2 = split_comma(elem) # split at commas - user defined function
			    store.append((body2)) # stores the new element

		# object to return - header only or body only
		if header_only == False:
			return store
		else:
			return header.split(',')
	
	else: 
		# for blank runs
		err = "{} was empty and was not added".format(code) # return an error
		return err




###########################################
######### Input(s) and Parameters #########
###########################################

# project codes
modern = np.genfromtxt('../data/prj_modern.csv', delimiter=',', dtype=str)
ancient = np.genfromtxt('../data/prj_ancient.csv', delimiter=',', dtype=str)

modern = modern[1:3]
ancient = ancient[1:3]

###########################################
############### Wraggling #################
###########################################

#########################################
### gather the data from NCBI entraze ###
	# store the header seperatly, the infotable gathered from NCBI forces a header each time. Other methods resulted in headers in the middle of the data.

# Modern
store_modern = [] # initialise storage 
header_modern = get_runinfo(modern[0], header_only=True) # header

# gather info for each of the codes
for i in modern:
	res = get_runinfo(i) 
	if type(res) != str: # if the code ran without errors
		store_modern = store_modern + get_runinfo(i)
	else:
		print(res) # print the error message if it occured


# Ancient
store_ancient = [] # initialise storage
header_ancient = get_runinfo(ancient[0], header_only=True) # header

# gather info for each of the codes
for i in ancient:	
	res = get_runinfo(i) 
	if type(res) != str: # if the code ran without errors
		store_ancient = store_ancient + get_runinfo(i)
	else:
		print(res) # print the error message if it occured


####################
### write to csv ###

# modern
with open('../data/info_modern.csv', 'w') as file:
	writer = csv.writer(file, delimiter=',')
	writer.writerow(i for i in header_modern) # add header
	for i in store_modern:
		writer.writerow(i)

# ancient
with open('../data/info_ancient.csv', 'w') as file:
	writer = csv.writer(file, delimiter=',')
	writer.writerow(i for i in header_ancient) # add header
	for i in store_ancient:
		writer.writerow(i)


