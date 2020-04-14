#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-01-14
# Last Modified: 2020-01-31



""" Read in project codes and expand adding from suplementary table
 	check to ensure there are no duplicates, etc."""

###########################################
################# Modules #################
###########################################

import csv
import numpy as np
from copy import deepcopy as dc

###########################################
######### Input(s) and Parameters #########
###########################################

####################
### project list ###

project_list = [] # initialise the list

with open('../data/raw_data/papers_projects.csv', 'r') as file:
	reader = csv.reader(file)
	for row in reader:
		project_list.append(row) # add row to list

##################################
### secondary list of projects ###

loc_orlando = '../data/raw_data/{}'.format(project_list[1][-3]) # location of orlando bioproject codes for Origins paper

project_orlando_all = [] # all data
project_orlando = [] # just the project codes

with open(loc_orlando, 'r') as file:
	reader = csv.reader(file)
	for row in reader:
		project_orlando_all.append(row) # all data
		project_orlando.append(row[3]) # project code

###########################################
############### Wraggling #################
###########################################

project_orlando = np.unique(project_orlando) # unique codes only
project_orlando = np.delete(project_orlando, -1) # strip header element

##############################################
### list of metadata and new project codes ###

list_to_append = []
tmp = dc(project_list[1]) # temp object deep copy - prevent overwrite of project_list

for i in project_orlando:
	tmp[0] = i # update the project code whilst maintainng other metadata
	list_to_append.append(tmp[:]) # append, slice to prevent shallow copy

####################################
### update the all projects list ###

project_list = project_list + list_to_append # append newly generated list
del project_list[1] # remove blank Orlando element that has now been populated

### add a unique reference 
for i in range(0, len(project_list)):
	if i == 0:
		project_list[i].insert(0, 'Source_Ref') # header
	else:
		project_list[i].insert(0, i) # numerical reference


with open('../data/cleaned_data/papers_projects_update.csv', 'w') as file:
	writer = csv.writer(file, delimiter=',')
	for i in project_list:
		writer.writerow(i)

############################
### Check for duplicates ###
	# they are present. not an issues, as mutiple sources may have been used across papers

all_project = [i[1] for i in project_list] # list of BioProject codes
all_project = all_project[1:] # remove header
all_project.sort()

#len(all_project) == len(np.unique(all_project)) # duplicates are present

####################################################
### Group project numbers by modern and ancient ####
	# 'both' has ancient and modern genomes associated with the code

# modern
modern = [i[1] for i in project_list if i[2] =='modern' or i[2] == 'both'] # modern project codes

modern = np.unique(modern)
np.savetxt('../data/cleaned_data/prj_modern.txt', modern, fmt='%s') # save to txt


# ancient
ancient = [i[1] for i in project_list if i[2] =='ancient' or i[2] == 'both'] # ancient project codes

ancient = np.unique(ancient)
np.savetxt('../data/cleaned_data/prj_ancient.txt', ancient, fmt='%s')




