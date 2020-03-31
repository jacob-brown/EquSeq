#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-03-31
# Last Modified: 2020-03-31



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








