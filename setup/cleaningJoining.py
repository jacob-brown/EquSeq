#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-01-27
# Last Modified: 2020-04-17



""" convert txt to csv and update info tables to include useful information in supplemenary tables """

###########################################
################# Modules #################
###########################################

import pandas as pd
import os
import csv
import numpy as np

###########################################
################# Functions ###############
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

# original downloaded txt infotable files
files = os.listdir('data/raw_data/infotables_original')
files.remove('.DS_Store')

# Orlando supplementary info
	# spans several project IDs 
df_S1_Orl = open_csv('data/raw_data/supplementary_data_from_studies/TableS1.final.csv')

# reference table for joining orlando 
ref_orlando = open_csv('data/raw_data/reference_orlando.csv')

###########################################
################# Wraggling ###############
###########################################

### remove projects that we don't need to join ###
ref_orlando = [i for i in ref_orlando if i[1] != '']
orlando_bio_projects = [i[0] for i in ref_orlando[1:]] # unique list of bioprojects

###################################################
### break string in orlando table to allow join ###

### store and delete header ###
head_S1 = df_S1_Orl[0] 
del df_S1_Orl[0]

### break the string ###
for i in range(0, len(df_S1_Orl)):
	tmp = df_S1_Orl[i] # store as tmp
	attach = tmp[0].split('_') # break at underscore
	# when underscore is in the sample name 
	if len(attach) > 3: # should only be a length of 3
		attach[2] = attach[2]+'_'+ attach[3] # concatinate and update
		attach = attach[0:3] # strip the last element
	df_S1_Orl[i] = tmp+attach # append and update the df

### update the header ### 
head_S1[3] = 'BioProject' # update header in prep for join
head_S1[2] = 'Breed_Info' # update breed as not to cause same name issues with joins
head_S1 = head_S1 + [1, 2, 'ID_approx'] # append new headers

### convert to pandas df and run through below: open, join, save ###
df_S1_Orl = pd.DataFrame(df_S1_Orl, columns = head_S1)


##########################################
### join suplementary tables and clean ###
	# also save as new csv, keeping originals

for i in files:
	#if i ==  'PRJEB9139.txt':
		#import ipdb
		#ipdb.set_trace()
	################################
	### open and save file paths ###
	path_open = 'data/raw_data/infotables_original/{}'.format(i) # open path
	no_ext = i.strip('.txt') # strip txt extension
	path_save = 'data/cleaned_data/infotables_update/{}.csv'.format(no_ext) # save path

	########################
	### open the txt file ###
		# files are comma sep but .txt
	df = pd.read_csv(path_open) 
	
	################################################################
	### conditions for if a supplementary table should be joined ###
		# some contain additional helpful information


	########################
	### Orlando Table S1 ###
	if no_ext in orlando_bio_projects and no_ext != 'SRP012260': # ensure it is in the ref table and not a project giving issues
		biop = [c[1] for c in ref_orlando if c[0] == no_ext] # bioprojects only

		df = pd.merge(df, df_S1_Orl, how='left', left_on=['BioProject', biop[0]], right_on=['BioProject','ID_approx'])
	



	###################################
	### PRJEB31613 & Source_Ref = 16 ###
		# '1-mmc1.csv' and PRJEB31613 - Source_Ref = 16
	if no_ext == 'PRJEB31613':

		sup16 = pd.read_csv('data/raw_data/supplementary_data_from_studies/1-mmc1.csv')

		# only the desired columns
		sup16 = sup16[['Sample name', 'Registration number', 'Period', 'Age (years ago) ',
		       'Site', 'latitude', 'longitude', 'Country', 'Species', 'Sex']]

		# left join - update df as using same name to save below
		df = pd.merge(df, sup16, how='left', left_on='Alias', right_on='Sample name')

	###################################
	### PRJEB10854 & Source_Ref = 6	###
		# mod: TableS1.1.csv and ancient: TableS1.2.csv 
	if no_ext == 'PRJEB10854':

		#df6 = pd.read_csv('data/infotables_update/PRJEB10854.csv')
		sup6_1 = pd.read_csv('data/raw_data/supplementary_data_from_studies/TableS1.1.csv')
		sup6_2 = pd.read_csv('data/raw_data/supplementary_data_from_studies/TableS1.2.csv')

		# trim to important data
		sup6_2 = sup6_2[['Horse ID', 'Site name and coordinates', 'Region of Origin', 'Age',]]
		sup6_1 = sup6_1[['Horse ID', 'Geographical coordinates', 'Region of Origin']]
		sup6_1['Age'] = '' # add empty vector for union
		sup6_1.columns = sup6_2.columns # update names for union
		sup6 = pd.concat([sup6_1, sup6_2],ignore_index=True) # union

		# left join - update df as using same name to save below
		df = pd.merge(df, sup6, how='left', left_on='Alias', right_on='Horse ID')

		#import ipdb
		#ipdb.set_trace()

	######################################
	### clean and subset unwanted data ###
	df = df[df['Assay Type'] == 'WGS'] # whole genomes only
	df = df[df['Organism'] != 'Gadus morhua'] # remove cod 

	#####################
	### save the data ###
	df.to_csv(path_save, index=False)



