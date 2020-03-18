#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-01-22
# Last Modified: 2020-03-18



"""  """

###########################################
################# Modules #################
###########################################

#import pandas as pd
import numpy as np
import os
import csv

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
				store[i][6] = 'ancient' # update
			elif store[i][4] == 'Yakutian horse':
				store[i][6] = 'modern' # update	
			store[i][7] = 'both_updated' # note the update


# end iteration # 

### re-add the header ###
store.insert(0, select_list)

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

# strip run and data_type in prep for groupnby
grp_by_vars = [i for i in head[1:] if i !='data_type'] 

# count number of runs and concatenate data_type
grp_df = df_store.groupby(grp_by_vars, as_index = False).agg({'data_type': ','.join, 'Run' : 'count'})

### only unique values in datatype ###
# pull out datatypes
df_data_type = grp_df.data_type.tolist()

# iterate through removing duplicates
for i in range(0,len(df_data_type)):
	string = df_data_type[i]
	tmp = string.split(',') # break into list
	uniq = set(tmp) # remove duplicates
	df_data_type[i] = ','.join(uniq) # rejoin and update df

# update pandas df
grp_df.data_type = df_data_type

# update column name
grp_df.rename(columns={'Run':'run_count'}, inplace=True)

# sort df
grp_df = grp_df.sort_values(by=['era', 'species', 'sub_group', 'data_type', 'BioProject', 'BioSample'], ascending=False)

# rearrange column order and remove age
grp_df = grp_df[['BioSample', 'BioProject', 'species', 'sub_group', 'era', 'data_type', 'run_count','comments']]

# save to csv
grp_df.to_csv('../data/cleaned_data/info_individual_grouped.csv', index=False)


############################################
### summarise at group and species level ###

# group by 
grp_pop_level = grp_df.groupby(['species', 'sub_group', 'era'], as_index = False).count()

# select only some vars
grp_pop_level = grp_pop_level[['species', 'sub_group', 'era', 'run_count']]

# rename count
grp_pop_level.rename(columns={'run_count':'individ_count'}, inplace=True)

# sort
grp_pop_level = grp_pop_level.sort_values(by=['era', 'species', 'sub_group'], ascending=False)

# save
grp_pop_level.to_csv('../data/cleaned_data/info_pop_grouped.csv', index=False)




