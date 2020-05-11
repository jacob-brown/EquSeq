#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-05-06
# Last Modified: 2020-05-07



"""  """

###########################################
################# Modules #################
###########################################

import os
import re 
from glob import glob
import numpy as np

###########################################
############## Function(s) ################
###########################################
def saveTxt(dirfile, listToSave, sep='\n'):
	with open(dirfile, 'w') as f:
		for num, val in enumerate(listToSave):
			if num == len(listToSave) - 1:
				f.write(val)
			else:
				f.write(val + sep)



def beagleMatch(inBeagle, snpList):
	# unzip - if required 
	if os.path.exists(inBeagle):

		os.system("gunzip {}".format(inBeagle))

	name = os.path.splitext(inBeagle)[0]
	basename = os.path.basename(name)

	f = open(name, 'r') 
	beagle = np.array(f.read().splitlines())

	### snp list ### 
	# allow passing of np.arrasy and lists 
	if type(snpList) == np.ndarray or type(snpList) == list:
		snps = snpList
	else:
		f_s = open(snpList, "r")
		snps = np.array(f_s.read().splitlines())

	header = beagle[0]

	### iterate through beagle codes and sample snps ###
	pattern = '[^\\t]*' # everything before first tab

	pos = lambda x: re.match(pattern, x).group().replace("_", ":")
	beagle_pos = np.array([pos(i) for i in beagle])

	snp_matches, b_in, s_in = np.intersect1d(beagle_pos, snps, \
								return_indices=True, assume_unique=True)

	# get the full beagle string 
	beagle_full = beagle[b_in]
	beagle_full = np.insert(beagle_full, 0, header)

	# update the snp list - removing matches 
	snps_no_match = np.delete(snps, s_in)

	print("\nsnps left to match: " + str(len(snps_no_match)) + "\n")

	### zip the original file ###
	os.system("gzip {}".format(name))

	# return match list in beagle format 
	return beagle_full, snps_no_match

###########################################
######### Input(s) and Parameters #########
###########################################

#inPath = "data/ancestry/test_success_beagle/*.beagle*"
inPath = "/rds/general/user/jb1919/ephemeral/ancestry/tmp_ALL_FILES/*.beagle*"
#outFile = "data/ancestry/test_success_beagle/ALL.modified.beagle"
outFile = "/rds/general/user/jb1919/ephemeral/ancestry/ancestry/ALL.modified.beagle"
#snplist = "data/ancestry/snp.list"
snplist = "/rds/general/user/jb1919/home/genomics/HorseGenomics/data/ancestry/snp.list"

### beagle files ###
beagle_ALL = glob(inPath)

# initiate 
#snp_count = sum(1 for line in open(snplist))
beagle_store = []

for file in beagle_ALL:
	beagles, snplist = beagleMatch(file, snplist)
	beagle_store.append(beagles)

	if len(snplist) == 0:
		print("matches complete, exiting.")
		break

	if len(beagles) == 1:
		print("no matches, moving to next file.")


### write ### 
saveTxt(outFile, beagle_store)


# add a better appedn np








