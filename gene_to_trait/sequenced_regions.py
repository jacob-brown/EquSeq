#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-04-22
# Last Modified: 2020-04-24



""" check bam file for regions sequenced, 
	comparing to desired variants """

###########################################
################# Modules #################
###########################################

import pysam
import mysql.connector
import re

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


###########################################
######### Input(s) and Parameters #########
###########################################

db = mysql.connector.connect(
	host="localhost",
	username="root",
	passwd="Bamboo91",
	database="OMIA")

cursor = db.cursor()

sql_File = "gene_to_trait/variants.sql"

bamfile = pysam.AlignmentFile("data/processed_sequences/new.rg.bam", "rb")


###########################################
############### Wraggling #################
###########################################

sql_Qry = open(sql_File, "r").read() # read string
cursor.execute(sql_Qry) # execute qry
results = cursor.fetchall()

header = cursor.column_names

i_chr = header.index('chromosome')
i_pos = header.index('g_or_m')
i_type = header.index('variant_type_name')
i_phene =  header.index('phenotype')

### quick filter of the variants ### -UPDATE LATER!!!!!!!! 
store = []
for elem in results:
	if(elem[i_type] == 'missense' \
		and elem[i_pos] != '' \
		and elem[i_chr] != ''):
		tmp_Pos = int(re.search('\d+', elem[i_pos]).group())
		tmp_chr = int(elem[i_chr])
		store.append([tmp_chr, tmp_Pos, elem[i_phene]])


# sort store list
store.sort(key=lambda x: (x[0], x[1]))

depth_store = []
for variant in store:
	chrom_pos = "chr"+str(variant[0])
	pos = variant[1]
	
	# open bam and determine depth
	depth = 0 # return value to 0 prior to running
	
	for read in bamfile.fetch(chrom_pos, pos, pos+1):
		depth = len([i for i in read.get_reference_positions() if i == pos])

	variant.insert(len(variant), depth)
	depth_store.append(variant)



store_best = [i for i in depth_store if i[len(depth_store[0])-1]  > 0]
store_notSeq = [i for i in depth_store if i[len(depth_store[0])-1]  == 0]

### write ### 
position_String = ["chr" + str(i[0]) + " "  + str(i[1]) for i in store]

#saveTxt("data/gene_variants/snp.list", position_String)

###########################################
############### Analysis ##################
###########################################


for variant in store:
	#try:
	chrPos = "chr" + str(variant[0])

	for pileupcolumn in bamfile.pileup(chrPos, variant[1], variant[1]+1):
		if(pileupcolumn.pos == variant[1]):
			seq = pileupcolumn.get_query_sequences(mark_matches=True, mark_ends=True, add_indels=True)
			mapQ = pileupcolumn.get_mapping_qualities()
			baseQ = pileupcolumn.get_query_qualities() # qualities
			depth = pileupcolumn.n
			print("\nposition: %s depth: %s seq: %s mapQ: %s base qual: %s " %
				(pileupcolumn.pos, str(depth), str(seq), str(mapQ), str(baseQ)))
	#except TypeError:
	#	pass




#bamfile.close()

###########################################
############### Plotting ##################
###########################################