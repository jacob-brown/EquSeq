#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-01-27
# Last Modified: 2020-01-28



""" convert txt to csv and update info tables to include useful information in supplemenary tables """

###########################################
################# Modules #################
###########################################

import pandas as pd
import os

###########################################
################# Functions ###############
###########################################

###########################################
######### Input(s) and Parameters #########
###########################################

# original downloaded txt infotable files
files = os.listdir('../data/infotables_original')
files.remove('.DS_Store')

###########################################
################# Wraggling ###############
###########################################


#########################################
### convert all files from txt to csv ###
for i in files:

	### open and save file paths ###
	path_open = '../data/infotables_original/{}'.format(i) # open path
	no_ext = i.strip('.txt') # strip txt extension
	path_save = '../data/infotables_update/{}.csv'.format(no_ext) # save path

	### open the txt file ###
		# files are comma sep but .txt
	df = pd.read_csv(path_open) 

	### conditions for if a supplementary table should be joined ###
		# some contain additional helpful information

	### PRJEB31613 & Source_Ref = 16 ###
		# '1-mmc1.csv' and PRJEB31613 - Source_Ref = 16
	if no_ext == 'PRJEB31613':

		#df16 = pd.read_csv('../data/infotables_update/PRJEB31613.csv')
		sup16 = pd.read_csv('../data/supplementary_data_from_studies/1-mmc1.csv')

		# only the desired columns
		sup16 = sup16[['Sample name', 'Registration number', 'Period', 'Age (years ago) ',
		       'Site', 'latitude', 'longitude', 'Country', 'Species', 'Sex']]

		# left join - update df as using same name to save below
		df = pd.merge(df, sup16, how='left', left_on='Alias', right_on='Sample name')

	### PRJEB10854 & Source_Ref = 6	###
		# mod: TableS1.1.csv and ancient: TableS1.2.csv 
	if no_ext == 'PRJEB10854':

		#df6 = pd.read_csv('../data/infotables_update/PRJEB10854.csv')
		sup6_1 = pd.read_csv('../data/supplementary_data_from_studies/TableS1.1.csv')
		sup6_2 = pd.read_csv('../data/supplementary_data_from_studies/TableS1.2.csv')

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
	### save the data ###
	df.to_csv(path_save)


