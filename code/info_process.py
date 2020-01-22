#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-01-21
# Last Modified: 2020-01-21



"""  """

###########################################
################# Modules #################
###########################################

import numpy as np
import csv

###########################################
############## Function(s) ################
###########################################

###########################################
######### Input(s) and Parameters #########
###########################################

### modern ####
modern = [] # initialise the list

with open('../data/manual_SraRunTable_modern.txt', 'r') as file:
	reader = csv.reader(file)
	for row in reader:
		modern.append(row) # add row to list


### ancient ####
ancient = [] 

with open('../data/manual_SraRunTable_ancient.txt', 'r') as file:
	reader = csv.reader(file)
	for row in reader:
		ancient.append(row) # add row to list


###########################################
############### Wraggling #################
###########################################

modern[0]
ancient[0]

###########################################
############### Analysis ##################
###########################################



###########################################
############### Plotting ##################
###########################################