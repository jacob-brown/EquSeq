#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-03-31
# Last Modified: 2020-04-06



""" Create a list of desired breeds,
	and pull the data. """

###########################################
################# Modules #################
###########################################

from fuzzywuzzy import fuzz, process
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


def write_csv(list_file, path):
	
	""" Write list to csv """

	with open(path, 'w') as f:
		writer = csv.writer(f, delimiter=',')
		for i in list_file:
			writer.writerow(i)

###########################################
######### Input(s) and Parameters #########
###########################################

# UK and Irish breeds
ukiebreeds = ['Shetland',
			'Eriskay',
			'Highland',
			'Clydesdale',
			'Fell ',
			'Dales',
			'Cleveland Bay',
			'English Thoroughbred',
			'Hackney',
			'Suffolk Punch',
			'Shire',
			'Welsh Mountain',
			'Welsh cob',
			'Welsh Pony',
			'Exmoor',
			'Dartmoor',
			'New Forest',
			'Connemara',
			'Irish cob',
			'Irish draught',
			'Irish sport',
			'Irish warmblood',
			'Kerry bog pony']


#individs = open_csv('../data/cleaned_data/info_individual_grouped.csv')
individs = open_csv('../data/cleaned_data/info_all.csv')

# store header
head_individs = individs[0]

###########################################
############### Wraggling #################
###########################################
# locations in header
grp_i = individs[0].index('sub_group')
era_i = individs[0].index('era')

# remove individuals with no breed data and modern only
individs = [i for i in individs if i[grp_i] !='' and i[era_i] == 'modern']

###########################################
############### Analysis ##################
###########################################

### store ratios for best matches of the string ###
store = []
for i in individs[1:]:
	ratios = process.extract(i[grp_i], ukiebreeds, limit=len(ukiebreeds)) 
	highest = process.extractOne(i[grp_i], ukiebreeds)
	tmp = [i[grp_i], highest[0], highest[1]]#, ratios]
	tmp = tmp + i
	
	# only store if the highest ratio is above 70
	if highest[1] > 70:
		store.append(tmp)



###########################################
############### Write csv #################
###########################################

# sort
store.sort(key = lambda x: x[2], reverse = True)

# add a header
header = ['sub_group', 'guess_grp', 'ratio'] + head_individs
store.insert(0, header)

# write to csv
write_csv(store, '../data/cleaned_data/breed_ratios.csv')


#----- save a list to use TEMP ------#

# headers to keep
list_head = ['sub_group', 'Run', 'BioSample', 'DATASTORE_filetype']

# list of index positions
list_i = [store[0].index(i) for i in list_head]

# extract the values
store_filter = []

for elem in store:
	tmp = [elem[i] for i in list_i]
	store_filter.append(tmp)


### get a couple of extra sampels as out groups ###
# use BioSample

####### remove later #############
out_samples = ['SAMEA104357346', 'SAMEA3504070']
out_info = [i for i in individs if i[1] in out_samples]

tmp_store = []
for elem in out_info:
	tmp = [elem[i] for i in [0, 1, 4, 5]]
	tmp = [tmp[2], tmp[0], tmp[1], tmp[3]]
	tmp_store.append(tmp)


store_filter.append(tmp_store[0])
store_filter.append(tmp_store[1])

# write to csv
write_csv(store_filter, '../data/cleaned_data/breed_sra_to_use.csv')

# write simple text file
runs = [i[1] for i in store_filter[1:]]
import numpy as np
np.savetxt('../data/cleaned_data/sra_runs.txt', runs, fmt='%s')
