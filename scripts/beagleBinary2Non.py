#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-10
# Last Modified: 2020-06-22



""" Convert binary beagle file to one with letters for alleles """

# python3 scripts beagleBinary2Non.py beaglein beagleout

###########################################
################# Modules #################
###########################################

import sys
import subprocess
import numpy as np

###########################################
############## Function(s) ################
###########################################

def saveTxt(dirfile, listToSave, sep='\n'):
	with open(dirfile, 'w') as f:
		for num, val in enumerate(listToSave):
			if num == len(listToSave) - 1:
				f.write(val)
			else:
				f.write(val + sep)

def beagleLetters(fileIn, fileOut):

	data = []
	task = subprocess.Popen(["cat", fileIn], stdout=subprocess.PIPE)
	for line in task.stdout:
		tmp = line.decode()
		tmp = tmp.replace("\n", "").split("\t")
		data.append(np.array(tmp))
	data = np.array(data)
	data

	bases = ["A", "C", "G", "T"]

	for val, elem in enumerate(data):
		if "allele" not in elem[1]:
			data[val][1] = bases[int(elem[1])]
			data[val][2] = bases[int(elem[2])]


	list_str = ["\t".join(i) for i in data]
	saveTxt(fileOut, list_str)


###########################################
################### Main ##################
###########################################



def main(argv):

	""" Main entry point run command """ 
	beagleLetters(argv[1], argv[2])
	return 0


if __name__ == "__main__":

    """Makes sure the "main" function is called from the command line"""
    
    status = main(sys.argv) 
    sys.exit(status)




