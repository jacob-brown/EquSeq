#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-01-22
# Last Modified: 2020-01-27



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
run_files = os.listdir('../data/infotables_update')
run_files.remove('.DS_Store')
run_files_no_ext = [i.strip('.csv') for i in run_files] # remove extension

### reference table ####
#ref_tab = pd.read_csv('../data/reference_infotable.csv')

ref_tab = open_csv('../data/reference_infotable.csv')

###########################################
############### Wraggling #################
###########################################

# loop through bioprojects, determine which column to gather and 
	# pull out that column
	# exclude any that are blank

#ref_tab[0] # headers
#ref_tab[0][3] # subgroup


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
	file_in = open_csv('../data/infotables_update/{}.csv'.format(pid)) 
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
	era_guess = ref_prj[5]
	comments = ref_prj[7]

	# iterate through each element (sequence run) of the filein
	for k in file_in:

		### pull data for each category - referencing location
		runs = k[0] # run ID
		sample = k[i_sample] # sample ID 
		species = k[i_species] # species 
		sub_group = k[i_sub_group] # how is the data grouped? population, breed, etc.
		data_type = k[i_data_type] # type of data - fasta, sam, etc.
		#age = k[i_age] # age

		#if pid == 'PRJEB10098':
		#	import ipdb
		#	ipdb.set_trace()

		# append to a new store list
			# iterate through each location and add the corresponding value to the list	
		tmp_store = [runs, sample, pid, species, sub_group, data_type, era, era_guess, '', comments]

		# append to master store
		store.append(tmp_store)


### remove Gadus morhua if present ###
	# present in one study
store = [i for i in store if i[3] != 'Gadus morhua'] 

### append the best guess of the age ###
for i in store:
	if i[6] != '':
		i.append(i[6])
	else:
		i.append(i[7])



# add the headers 
head = ['Run', 'BioSample', 'BioProject', 'species', 'sub_group', 'data_type', 'era', 'era_guess', 'age', 'comments', 'era_best_guess']


df_store = pd.DataFrame(store, columns=head)


df_store.to_csv('../data/run_info_summary.csv')




### rename to csv's ###
	# quickly viewing become easier

if False:

	for i in files:
		og = '../data/infotables_update/{}'.format(i)
		strp = i.strip('.txt') + '.csv'
		rnam = '../data/infotables_update/{}'.format(strp)
		os.rename(og, rnam) 




	df_modern.columns

	### remove non Equus species ###
	df_ancient = df_ancient[df_ancient.Organism != 'Gadus morhua']

	### remove RNA-Seq ###
	df_modern = df_modern[df_modern['Assay Type'] != 'RNA-Seq']

	### check for duplicates ###

	duplicates = df_modern[df_modern.Run.isin(df_ancient.Run)]

	# which project IDs are the duplicates found in 
	dup_bioPrj = duplicates.BioProject.unique() # bio project
	dup_run = duplicates.Run.unique() # run ids

	# print duplication project codes
	if len(duplicates) != 0:
		print("Duplicates found in {}, check the data.".format(dup_bioPrj))


	########
	### Resolve duplication
		# much of the duplication in bioproject codes is due to 
		# grouping of ancient and modern in the same sets

	# import data for splitting  
	df_modern_split = pd.read_csv('../data/supplementary_data_from_studies/TableS1.1.csv')
	df_ancient_split = pd.read_csv('../data/supplementary_data_from_studies/TableS1.2.csv')


	'PRJEB10854' 
	'PRJEB19970'







