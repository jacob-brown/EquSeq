#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-01-22
# Last Modified: 2020-03-31



""" Summarise inforun tables into:
	1) all data
	2) individual sample level
	3) population/breed and age level """

###########################################
################# Modules #################
###########################################

#import pandas as pd
import numpy as np
import os
import csv
import pandas as pd

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


def duplicate_runs(data):

	""" check for duplicate runs """

	runs = [i[0] for i in data]

	if len(np.unique(runs)) != len(data): # should be only uniques 
		return 'Duplicate runs found.'
	else:
		return 0


def no_group(data):
	
	""" Group/Breed info not present.
		Return ProjectIDs """

	tmp = [i[2] for i in data if i[4]=='']
	return np.unique(tmp)

def issue_BioProjects(all_files, ID_keep, ID_ignore, \
						out='../data/cleaned_data/issue_BioProject.txt'):
	
	""" save file of BioProjectIDs that still need work.
		those that are present but are not used """

	# strip 
	files_no_ext = [i.strip('.csv') for i in all_files]

	### check that all projects are present, unless noted to ignore ###
		# save a file of those that still need work

	issue_IDs = [] # tmp storage

	for i in files_no_ext:
		if i not in ID_keep: # project ID not in project list
			issue_IDs.append(i)
			if i not in ID_ignore: # unless specified
				print('{}.csv present, but not specified in reference table'.format(i))

	# save the IDs that need further exploration
	np.savetxt(out, issue_IDs, delimiter='', fmt='%s')

	return 'Issue BioProject IDs saved.'

###########################################
######### Input(s) and Parameters #########
###########################################


### runinfo table files ###
run_files = os.listdir('../data/cleaned_data/infotables_update')


# remove DS_Store if present
if '.DS_Store' in run_files:
	run_files.remove('.DS_Store')

# reference table - details columns of interest 
ref_tab_all = open_csv('../data/raw_data/reference_infotable.csv')

# remove the ignore=yes projects and strip ignore column 
ref_tab = [i[0:6] for i in ref_tab_all if i[6] != 'y'] 


### project IDs ### 
prj_ID = [i[0] for i in ref_tab] # Project IDs to be used
del prj_ID[0] # remove header
prj_ID_ignored = [i[0] for i in ref_tab_all if i[6] == 'y'] # IDs ignored


###########################################
############### Wraggling #################
###########################################


### filter columns of infotables ###
	# allowing for unions and summarising
store = []

for pid in prj_ID:

	# open file
	file_in = open_csv('../data/cleaned_data/infotables_update/{}.csv'\
		.format(pid)) 
	header = file_in[0] # header list
	del file_in[0] # strip header

	### present in all
	# SRA Study in all, not always SRA_ass
	# DATASTORE filetype in all
	# Organism in all

	### columns of interest, from reference table ###
		# ['BioProject', 'sub_group', 'era', 'era_guess', 'age', 'comments']
	ref_pivot = [i for i in ref_tab if i[0] == pid][0] # remove nesting with [0] 

	### Extract information of interest and append to data ###
		

	# sub_group - find the index within the header
	i_sub = header.index(ref_pivot[1]) 

	# append (and extend, when list) subgroup and reftab info
	for elm in file_in:
		elm.append(elm[i_sub]) # sub group
		elm.extend(ref_pivot[2:6]) # era, era_guess, age, comments


	#---- write to function ----#
	# update header	- unsure same order as appended above
	header.extend(['sub_group', 'era', 'era_guess', 'age', 'comments'])

	### columns to keep ###
	select_list = ['Run', 'BioSample', 'BioProject','Organism',\
					'sub_group', 'DATASTORE filetype', 'SRA Study',\
					'era', 'era_guess', 'age', 'comments']

	# find their position within the header
	index_list = [header.index(i) for i in select_list]
	
	# pull from the main data
	for row_in in file_in:
		tmp_selected = [row_in[i] for i in index_list] 
		store.append(tmp_selected)


	### update 'both' field for Yakutian horse - PRJEB10854 ###
		# hardcode for finer detail
	for i in range(0, len(store)):
		if store[i][7] == 'both': # check both
			if store[i][4] == 'Ancient horse from Yakutia': 
				store[i][7] = 'ancient' # update
			elif store[i][4] == 'Yakutian horse':
				store[i][7] = 'modern' # update	


# end iteration # 

### re-add the header ###
select_list = [i.replace(' ','_') for i in select_list] # replace spaces
store.insert(0, select_list) # header to store


### save to csv ###
with open('../data/cleaned_data/info_all.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(store)

##############
### Checks ###

# 0 notes no duplicates
duplicate_runs(store) 

# blank entries
no_group(store) 

# bioproject IDs that need attention - saves file
issue_BioProjects(run_files, prj_ID, prj_ID_ignored)

##################################################
### summarise the data into individual samples ###

# convert to pandas df for easy group by
store_df = pd.DataFrame(store[1:],columns=store[0])


#--- individual sample grouping ---# 
	# group by everything except Run and DATASTORE_filetype
	# concatenate DATASTORE_filetype

# variables to group by
group_vars_ind = np.delete(store[0], [0, 5]) # head except run and DATASTORE_filetype
group_vars_ind = list(group_vars_ind) # ensure non-numpy
# group by
ind_df = store_df.groupby(group_vars_ind, as_index = False)\
			.agg({'DATASTORE_filetype': ','.join, 'Run' : 'count'})

# correct the datastore file type
dt_array = []
for i in np.array(ind_df.DATASTORE_filetype):
	tmp = i.split(',')
	tmp = np.unique(tmp)
	tmp = str(tmp).strip("[ ]")
	dt_array.append(tmp)

# update main variable
ind_df.DATASTORE_filetype = dt_array

# update column name
ind_df.rename(columns={'Run':'run_count'}, inplace=True)

# save to csv
ind_df.to_csv('../data/cleaned_data/info_individual_grouped.csv', index=False)

#--- population grouping ---#
# variables to group by

pop_df = ind_df.groupby(['Organism', 'sub_group', 'era', 'age'], \
	as_index = False).count()

# select only some vars
pop_df = pop_df[['Organism', 'sub_group', 'era', 'age', 'run_count']]

# rename count
pop_df.rename(columns={'run_count':'individ_count'}, inplace=True)

# sort
pop_df = pop_df.sort_values(by=['era', 'Organism', 'sub_group'], \
	ascending=False)

# save
pop_df.to_csv('../data/cleaned_data/info_pop_grouped.csv', index=False)
