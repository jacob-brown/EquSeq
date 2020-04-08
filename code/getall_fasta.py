#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-04-07
# Last Modified: 2020-04-08

""" get raw fasta data from ENA """

###########################################
################# Modules #################
###########################################

import requests
import os 
import sys
import csv
#from optparse import OptionParser 

###########################################
############## Function(s) ################
###########################################

def getTableInfo(code, getall=True):

	""" get run info data from ENA """

	print("requesting infotable")

	# get all the run info data, or just the fastq_ftp info
	if getall: 
		link = "https://www.ebi.ac.uk/ena/data/warehouse/"\
		"filereport?result=read_run&accession={}".format(code)
	else:
		link = "https://www.ebi.ac.uk/ena/data/warehouse/"\
		"filereport?result=read_run&fields=fastq_ftp&accession={}".format(code)
	r = requests.get(link)
	text = r.text # single string text of table
	head_body = text.split('\n') # split into head and body
	# split into header and body, excluding redundent blank elements
	list_out = [i.split('\t') for i in head_body if i != '']
	return list_out

def getFasta(code, path_out):

	print("requesting fasta: " + code)

	info_run = getTableInfo(code, getall=False)

	# avoid blank infotable returns
	if len(info_run) > 1 :

		### construct the queries ###
		queries = info_run[1][0].split(';')

		for n, elem in enumerate(queries):

			base = os.path.basename(elem) # strip path
			# create the string
			#qry = "wget -O {}/{} {}".format(path_out, base, elem)
			qry = "curl -L {} -o {}/{}".format(elem, path_out, base)
			
			queries[n] = qry # update queries with new command


		### concatenate into single command and run ###
		command = ';'.join(queries)
		os.system(command) # run command
		
	else: 
		return 'Info table empty for: ' + code + '\n check code.'

def accessionCode(jobN, file):

	""" get the asscession code from the list """
	
	with open(file, 'r') as f:
		tmp = f.readlines() # read file
	
	code = tmp[jobN].strip('\n') # filter and strip new line
	
	return code

###########################################
################# main ####################
###########################################

def main(argv):

	""" Main entry point run command """ 
	
	### catch args ##
	job_num = int(os.environ['PBS_ARRAY_INDEX'])
	sra_list_loc = argv[1]	
	fasta_path_out = argv[2]
	
	# get accession code and run fasta retrival
	code = accessionCode(job_num, sra_list_loc)
	getFasta(code, fasta_path_out) # code and dirictory
	
	return 0


if __name__ == "__main__":

    """Makes sure the "main" function is called from the command line"""
    
    status = main(sys.argv) 
    sys.exit(status)

