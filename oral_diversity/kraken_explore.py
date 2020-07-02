#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-07-01
# Last Modified: 2020-07-02



""" count reads per taxaa based on a filtering threshold
	generates 2 csvs (pass and fail) """

# usage example
# python3 oral_diversity/kraken_explore.py -i sandbox/kraken/all.tmp.kraken -o sandbox/kraken/stats -t 0.2

###########################################
################# Modules #################
###########################################

import numpy as np
from itertools import compress 
import argparse
import csv
import subprocess

###########################################
############## Function(s) ################
###########################################

def write_csv(list_file, path):
	
	""" Write list to csv """

	with open(path, 'w') as f:
		writer = csv.writer(f, delimiter=',')
		for i in list_file:
			if isinstance(i, list):
				writer.writerow(i) # multi-column
			else:
				writer.writerow([i]) # single column
	
# sub process wrapper
def subProWrap(command, returnList=True):
	""" wrapper for a subprocess command 
		returns list as default """
	p = subprocess.Popen([command], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	out, er = p.communicate()
	if(returnList):
		out_string = out.decode()
		files = out_string.replace("\t", ":").split("\n")
		files.remove("")
		return files

def cqScore(elem):
	""" calculate kraken qc score for 
		read classification string """

	delem = [np.array(i.split(":")).astype(int) for i in elem]
	ids = [i[0] for i in delem]
	values = np.array([i[1] for i in delem])
	groups = np.unique(ids)
	sum_values = np.array([values[ids==i].sum() for i in groups])
	
	# check if 0 is present 
	if 0 not in groups:
		return 1 # all kmers have mapped to a taxa
	else:
		# which values are category 0 
		q = sum_values[groups == 0][0]
		# remove 0 fields
		qindex = list(sum_values).index(q)
		sum_values = np.delete(sum_values, qindex)
		#sum_values = sum_values[sum_values != q]
		groups = groups[groups != 0]
		
		# and which is the taxa with the highest match
		c = sum_values[list(sum_values).index(sum_values.max())]
		
		return c/q 

# Print iterations progress
def printProgressBar(iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()

# pass kraken file, calculate 
def krakenPass(kraken, out, thresh = 0.2):
	
	""" generate summary counts per taxa for kraken
		classification above a C/Q threshold (0-1) """

	print("reading kraken file")
	#counter_unclass = 0
	dict_pass = {}
	dict_fail = {"U" : 0} 

	print("getting line count")
	counterLine = 0
	totalLines = int(subProWrap("wc < "+ kraken +" | awk {'print $1'}")[0])
	print(str(totalLines) + " lines to read")
	print("\ncomputing cq scores and counting...\n")
	with open(kraken, "r") as infile:
		for line in infile:
			
			# print progress bar
			#printProgressBar(iteration=counterLine, total=totalLines, prefix = 'Progress:', suffix = 'Complete', length = 50)
			if counterLine % 50000 == 0:
				print(str(counterLine) + "/" + str(totalLines))

			counterLine = counterLine + 1
			# count unclassified
			if line[0] == "U":
				dict_fail["U"] += 1
			else:
				tmp_line = line.replace("\n", "").split("\t")
				split_class = tmp_line[4].replace("|:| ", "").split(" ")

				# passes threshold value
				
				if cqScore(split_class) >= thresh:

					if tmp_line[2] not in dict_pass.keys():
						# start the counter if a new taxa
						dict_pass[tmp_line[2]] = 1 
					else:
						# add a value to the taxa
						dict_pass[tmp_line[2]] += 1
				else:
					# lower than threshold 
					if tmp_line[2] not in dict_fail.keys():
						# start the counter if a new taxa
						dict_fail[tmp_line[2]] = 1 
					else:
						# add a value to the taxa
						dict_fail[tmp_line[2]] += 1



	print("\n"+str(sum(dict_pass.values())) + "/" + str(sum(dict_pass.values()) + sum(dict_fail.values())) + " reads remain\n")

	### save
	# pass 
	pass_fileout = out+".pass.csv"
	print("saving: " + pass_fileout)

	with open(pass_fileout, 'w') as csv_file:  
	    writer = csv.writer(csv_file)
	    for key, value in dict_pass.items():
	       writer.writerow([key, value])

	# failed
	fail_fileout = out+".fail.csv"
	print("saving: " + fail_fileout)

	with open(fail_fileout, 'w') as csv_file:  
	    writer = csv.writer(csv_file)
	    for key, value in dict_fail.items():
	       writer.writerow([key, value])




###########################################
################# Options #################
###########################################

parser = argparse.ArgumentParser(description=\
		'Generate new filtered kraken files using a confidence threshold')



# input file
parser.add_argument("-i", "--input", dest="infile", type=str,
					required=True, help="kraken input file")

#  output prefix
parser.add_argument("-o", "--out", dest="outprefix", type=str,
                  required=False, default="./", help="prefix of output files")

# threshold
parser.add_argument("-t", "--thresh", dest="thold", type=float,
                  required=True, default=0.2, \
                  help="threshold to filter c/q float from (0-1)")

# define args
args = parser.parse_args()

###########################################
################# main ####################
###########################################

# python3 oral_diversity/kraken_explore.py -i sandbox/kraken/tmp.kraken -o sandbox/kraken/pass -t 0.2

krakenPass(kraken=args.infile, out=args.outprefix, thresh=args.thold)



#fopen = open(kraken, "r")
#file = fopen.readlines()
#data = [i.replace("\n", "").split("\t") for i in file]
#split_class = [i[4].replace("|:| ", "").split(" ") for i in data]
#
## run functions and filter on cq over threshold
#
#for i in split_class:
#	tmp = cqScore(i)
#
#cqValues = [cqScore(i) for i in split_class]
#keepCond = [i >= thresh for i in cqValues] # boolean
#passed_reads = list(compress(file, keepCond)) # filter data
#
#print(str(len(passed_reads)) + "/" + str(len(data)) + " reads remain")
#
## save
#fileout = out+".kraken"
#print("saving: " + fileout)
#saveTxt(fileout, passed_reads, sep='')



## C/Q
#	# C is the number of k-mers mapped to LCA values in the clade rooted at the label
#	# Q is the number of k-mers in the sequence that lack an ambiguous nucleotide
#
#### single taxa ###
#'0:28 109376:3 0:35 |:| 109376:41 0:25'
##per taxa: classified/unclassified 
#
#(3+41)/(28+35+25) # classification is 50% confident, use 20%
#
#### multi-taxa ###
#'0:41 9606:2 0:23 |:| 0:8 9606:5 0:1 131567:5 9606:5 131567:7 1:5 131567:3 1:2 131567:1 1283:1 9606:4 0:19'
#q=(41+23+8+1+19)
## 9606
#(2+5+5+4)/q
## 131567
#(5+7+3+1)/q
## 1283
#(1)/q
#
## this would be removed with 20% confidence