#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-01-22
# Last Modified: 2020-03-09



"""  """

###########################################
################# Modules #################
###########################################

import pandas as pd
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

###########################################
######### Input(s) and Parameters #########
###########################################


### runinfo table files ###
run_files = os.listdir('../data/cleaned_data/infotables_update')

# remove DS_Store if present
if '.DS_Store' in run_files:
	run_files.remove('.DS_Store')

# strip csv extension
run_files_no_ext = np.array([i.strip('.csv') for i in run_files]) 

### reference table ####
#ref_tab = pd.read_csv('../data/reference_infotable.csv')

ref_tab = open_csv('../data/raw_data/reference_infotable.csv')
ref_tab = [i[0:8] for i in ref_tab if i[8] != 'y'] # remore the ignore projects and strip ignore column 

###########################################
############### Wraggling #################
###########################################

# loop through bioprojects, determine which column to gather and 
	# pull out that column
	# exclude any that are blank

### seperate groups with and without information ###
all_PrjID = [i[0] for i in ref_tab] # all project ids, keep header as using for index
no_df_groups = [i for i in ref_tab if i[3] == ''] # without sub_group
df_groups = [i for i in ref_tab if i[3] != ''] # with sub_group
prj_ID = [i[0] for i in df_groups] # Project IDs to be used
del prj_ID[0] # remove header

### use the reference table to pull out desired columns and append
store = []

for pid in prj_ID:

	# open file
	file_in = open_csv('../data/cleaned_data/infotables_update/{}.csv'.format(pid)) 
	header = file_in[0] # header list
	del file_in[0] # strip header

	# find location of project code in the ref table
	prj_loc = all_PrjID.index(pid)
	ref_prj = ref_tab[prj_loc] # refernece for the project of interest

	# find location with the header of the read in file
	i_sample = header.index('BioSample')
	i_sub_group = header.index(ref_prj[3]) # sub_group index (breed, location, etc.)
	i_data_type = header.index(ref_prj[1]) # data type (fasta, bam, sam, etc.)
	i_species = header.index(ref_prj[2]) # species 
	#i_age = header.index(ref_prj[6]) # age 

	### ref table info ###
	era = ref_prj[4]
	#era_guess = ref_prj[5]
	comments = ref_prj[7]

	# iterate through each element (sequence run) of the filein
	for k in file_in:

		### pull data for each category - referencing location
		runs = k[0] # run ID
		sample = k[i_sample] # sample ID 
		species = k[i_species] # species 
		sub_group = k[i_sub_group] # how is the data grouped? population, breed, etc.
		data_type = k[i_data_type] # type of data - fasta, sam, etc.

		# append to a new store list
			# iterate through each location and add the corresponding value to the list	
		tmp_store = [runs, sample, pid, species, sub_group, data_type, era, '', comments]

		# append to master store
		store.append(tmp_store)




### note the BioProjects with missing data ###
#issue_IDs = np.unique([i[2] for i in store if i[4] == '']) # those that had partial joins
#passed_IDs = np.unique([i[2] for i in store]) # ids that passed
#add_issues = run_files_no_ext[np.isin(run_files_no_ext, passed_IDs, invert=True)] # list those values not in the main df i.e. no info
#issue_IDs = np.append(issue_IDs, add_issues) # append partial and failed ids


### tidy up the df prior to save ###
store = [i for i in store if i[4] != ''] # no blank entries

### update 'both' field for Yakutian horse - PRJEB10854 ###
	# hardcode for finer detail
for i in range(0, len(store)):
	if store[i][6] == 'both': # check both
		if store[i][4] == 'Ancient horse from Yakutia': 
			store[i][6] = 'ancient' # update
		elif store[i][4] == 'Yakutian horse':
			store[i][6] = 'modern' # update

# add the headers 
head = ['Run', 'BioSample', 'BioProject', 'species', 'sub_group', 'data_type', 'era', 'age', 'comments']

df_store = pd.DataFrame(store, columns=head)
df_store.to_csv('../data/cleaned_data/info_all.csv', index=False)

### which bioprojects failed or still need work ###
# unique list of projectIDs
proj_all = open_csv('../data/cleaned_data/papers_projects_update.csv')
proj_all = [i[1] for i in proj_all] # projectIDs only
del proj_all[0] # strip header (not an ID)
proj_all = np.unique(proj_all) # unique only

# IDs that have been processed
passed_ID = np.unique(df_store.BioProject)

# compare the processed with all
issue_IDs = proj_all[np.isin(proj_all, passed_ID, invert=True)] 
issue_IDs = issue_IDs.astype(str)
np.savetxt('../data/cleaned_data/issue_BioProject.txt', issue_IDs, delimiter='', fmt='%s')


### check for duplicates ###
if len(np.unique(df_store.Run)) != len(df_store):
	print('Duplicate runs found.')

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




