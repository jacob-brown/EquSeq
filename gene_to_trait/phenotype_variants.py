#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-04-22
# Last Modified: 2020-06-12

""" regions to snp call. based on omia and QTLdb """

###########################################
################# Modules #################
###########################################

import mysql.connector
import re
import hgvs.parser
import csv
import numpy as np
import getpass

###########################################
############## Function(s) ################
###########################################

### save a text file without a new line at the end
def saveTxt(dirfile, listToSave, sep='\n'):
	with open(dirfile, 'w') as f:
		for num, val in enumerate(listToSave):
			if num == len(listToSave) - 1:
				f.write(val)
			else:
				f.write(val + sep)

### return the important info from hgvs class ###

def hgvsGetInfo(SequenceVariant):
	# start, end, type, ref, alt, 
	return [SequenceVariant.posedit.pos.start.base,
				SequenceVariant.posedit.pos.end.base, \
				SequenceVariant.posedit.edit.type, \
				SequenceVariant.posedit.edit.ref, \
				SequenceVariant.posedit.edit.alt]

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

print("sql password...")
db = mysql.connector.connect(
	host="localhost",
	username="root",
	passwd=getpass.getpass(),
	database="OMIA")

cursor = db.cursor()
sql_File = "gene_to_trait/variants.sql"
omia_catch = open_csv('data/gene_variants/omia_catch_all.csv')

###########################################
############### Wraggling #################
###########################################
print("accessing sql")

sql_Qry = open(sql_File, "r").read() # read string
cursor.execute(sql_Qry) # execute qry
results = cursor.fetchall()

header = cursor.column_names

### update header names ###
match_list = ["g_or_m", "c_or_n", "p"]
replace_list = ["coordinate", "coding_position", "protein_position"]
header = list(header) 

for i in range(0, len(match_list)):
	i_loc = header.index(match_list[i])
	header[i_loc] = replace_list[i]

i_chr = header.index('chromosome')
i_coord = header.index('coordinate')
i_phene =  header.index('phenotype')
i_varID = header.index('variant_id')
i_breed = header.index('breed_name')



### decode hgvs ###
hp = hgvs.parser.Parser()

# some genomic (.g) coordinates are not correctly formatted or are empty
	# as cDNA and/or protein positions are given (.c .p)
store = []
for var in results:
	if var[i_coord] != ''and var[i_chr] != '':
		try:
			# parser requires an accession number
			tmp_coor =  "accN:" + var[i_coord] 
			hgvs_var = hp.parse_hgvs_variant(tmp_coor)
			info =  hgvsGetInfo(hgvs_var)
			
			to_append = [int(var[i_varID]), int(var[i_chr]), var[i_phene], var[i_breed]]
			to_append.extend(info)
			store.append(to_append)
		except:
			pass


# varID, chr, phen, breed, position, position_stop, type, reference_nucleotide, new_nucleotide 

### filter omio csv - SQL db dump was not updated on OMIAs side ###
	# hgvs code is the most unique for the presence/absence check
hgvs_codes = [i[8] for i in results if i[8] != '']

### run parser on omia caught data ###

store_omia_catch = []

for o in omia_catch:
	if o[9] != '' and o[10] != '' \
	and o[10].replace(',','') not in hgvs_codes:
		tmp = o
		tmp[10] = tmp[10].replace(',','')
		
		try:
			# parser requires an accession number
			tmp_coor =  "accN:" + var_o[10] 
			hgvs_var = hp.parse_hgvs_variant(tmp_coor)
			info =  hgvsGetInfo(hgvs_var)
			
			to_append = [var_o[0], int(var_o[9]), var_o[3], var_o[2]]
			to_append.extend(info)
			store_omia_catch.append(to_append)
		except:
			pass


### append missing omia sites ###
store = store + store_omia_catch

# sort store list - chr and start position
store.sort(key=lambda x: (x[1], x[4]))

### write ### 
position_String = ["chr" + str(i[1]) + "\t"  + str(i[4]) + "\t"  + str(i[5]) for i in store]

# angsd list format chr:pos
print("making position list files")

store_str = []
for i in store:
	for c in range(i[4]-1, i[5]+1):
		tmp = str("chr" + str(i[1]) + ":" + str(c))
		store_str.append(tmp)

split_str = np.array_split(np.array(store_str), 10)

for val, elem in enumerate(split_str):
	saveTxt("data/gene_variants/trait.snps/trait.snp.{}.list".format(val), elem)


# bed list format chr start end - samtools usage
print("making bed files")

store_bed = []
for i in store:
	tmp = str("chr" + str(i[1]) + "\t" + str(i[4]-1) + "\t" + str(i[5]+1))
	store_bed.append(tmp)

split_bed = np.array_split(np.array(store_bed), 10)

for val, elem in enumerate(split_bed):
	saveTxt("data/gene_variants/trait.snps/trait.snp.{}.bed".format(val), elem)


# non bed list
tab_sep = [i.replace(":", "\t") for i in store_str]
split_tab = np.array_split(np.array(tab_sep), 10)

for val, elem in enumerate(split_tab):
	saveTxt("data/gene_variants/trait.snps/trait.snp.{}.tab".format(val), elem)

# write general info csv
write_csv(store, "data/gene_variants/snp.list.info.csv")










