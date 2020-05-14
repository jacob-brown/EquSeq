#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-04-14
# Last Modified: 2020-05-11



""" helpful functions defined by Jacob Brown """


# timer
import time
start = time.time() # start the timer from import

def timer():
	
	end = time.time()
	duration = end-start
	duration = round(duration)
	string = "\n..........................\n"\
			"   Time elapsed: {} sec"\
			"\n..........................\n"\
			.format(duration)

	print(string)


### save a text file without a new line at the end
def saveTxt(dirfile, listToSave, sep='\n'):
	with open(dirfile, 'w') as f:
		for num, val in enumerate(listToSave):
			if num == len(listToSave) - 1:
				f.write(val)
			else:
				f.write(val + sep)



def open_csv(file):
	import csv
	""" open a csv into a list format """

	tmp = [] # initialise the list
	with open(file, 'r') as f:
		reader = csv.reader(f)
		for row in reader:
			tmp.append(row) # add row to list

	return tmp