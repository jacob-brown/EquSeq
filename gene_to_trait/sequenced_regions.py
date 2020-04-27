#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-04-22
# Last Modified: 2020-04-27



""" check bam file for regions sequenced, 
	comparing to desired variants """

###########################################
################# Modules #################
###########################################

import pysam
import mysql.connector
import re
import hgvs.parser

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
			
			to_append = [int(var[i_varID]), int(var[i_chr]), var[i_phene]]
			to_append.extend(info)
			store.append(to_append)
		except:
			pass


# sort store list - chr and start position
store.sort(key=lambda x: (x[1], x[3]))


depth_store = []
for variant in store:
	chrPos = "chr"+str(variant[1])
	
	# open bam and determine depth
	depth = 0 # return value to 0 prior to running
	
	for read in bamfile.fetch(chrPos, variant[3], variant[4]+1):
		depth = len([i for i in read.get_reference_positions() if i >= variant[3] and i <= variant[4] ])

	variant.insert(len(variant), depth)
	depth_store.append(variant)


store_best = [i for i in depth_store if i[10] > 0]
store_notSeq = [i for i in depth_store if i[10]  == 0]

### write ### 
#position_String = ["chr" + str(i[0]) + " "  + str(i[1]) for i in store]

#saveTxt("data/gene_variants/snp.list", position_String)

###########################################
############### Analysis ##################
###########################################
variant = [i for i in store if i[0] == 890][0]

for variant in store:
	
	chrPos = "chr" + str(variant[1])
	
	for pileupcolumn in bamfile.pileup(chrPos, variant[3], variant[4]+1):
		if(pileupcolumn.pos >= variant[3] \
			and pileupcolumn.pos <= variant[4]):
			#for pileupread in pileupcolumn.pileups:
			#	tmp_info = [pileupread.is_del, \
			#					pileupread.query_position, \
			#	pileupread.alignment.query_sequence[pileupread.query_position]]
			#					
			#	print(tmp_info)
			seq = pileupcolumn.get_query_sequences(mark_matches=True, \
				mark_ends=True, add_indels=True)
			mapQ = pileupcolumn.get_mapping_qualities()
			baseQ = pileupcolumn.get_query_qualities() # qualities
			depth = pileupcolumn.n

			print("\nposition: %s depth: %s seq: %s mapQ: %s base qual: %s " %
				(pileupcolumn.pos, str(depth), str(seq), str(mapQ), str(baseQ)))
			print(variant[2] + ", "+variant[5])
			print("ref: " + variant[6] + " " + str(seq))
			print("var: " + str(variant[7]) + " " + str(seq))





#for pileupcolumn in bamfile.pileup(chrPos, variant[3], variant[4]+1):
#    print ("\ncoverage at base %s = %s" %
#           (pileupcolumn.pos, pileupcolumn.n))
#    for pileupread in pileupcolumn.pileups:
#        if not pileupread.is_del and not pileupread.is_refskip:
#            # query position is None if is_del or is_refskip is set.
#            print ('\tbase in read %s = %s' %
#                  (pileupread.alignment.query_name,
#                   pileupread.alignment.query_sequence[pileupread.query_position]))
#


#bamfile.close()

###########################################
############### Plotting ##################
###########################################