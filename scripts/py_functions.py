#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-04-14
# Last Modified: 2020-06-26



""" helpful functions defined by Jacob Brown """

# timer
import time

def timer(state="S"):
	global startTimer

	"""" timer(), timer("e") to start and stop timer """
	if state.lower()=='s':
		
		print("timer started")
		startTimer = time.time() # start the timer from import	
	
	elif state.lower()=='e':

		duration = time.time() - startTimer
		duration = round(duration)
		string = "\n..........................\n"\
				"   Time elapsed: {} sec"\
				"\n..........................\n"\
				.format(duration)

		print(string)
	
	else:
		stop("state not found, s or e only.")


# write csv
def write_csv(list_file, path):
	
	""" Write list to csv """

	with open(path, 'w') as f:
		writer = csv.writer(f, delimiter=',')
		for i in list_file:
			if isinstance(i, list):
				writer.writerow(i) # multi-column
			else:
				writer.writerow([i]) # single column
	


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
