#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-05-04
# Last Modified: 2020-05-06



""" Generate sites for use in 
	genotype liklihood and PCA  """

###########################################
################# Modules #################
###########################################

import pysam
import pickle
import numpy as np
import time
from natsort import natsorted
import argparse
import pandas as pd 


###########################################
################# Options #################
###########################################

parser = argparse.ArgumentParser(description=\
		'Generate snp position list on a fuzzy position of a bam file.')


# command 
choose = ["snps", "subsample", "split"]
parser.add_argument("-c", "--command", dest="command", type=str,
					required=True, choices=choose, 
					help="choice function. snps: generate pickle of all sites; subsample: "\
						"the sites; chrmSplit: split outfiles by chromosome.")

# input bam file
parser.add_argument("-i", "--input", dest="infile", type=str,
					required=True, help="bam IN_BAM to have snps sampled from", \
					metavar="IN_BAM")

#  output directory
parser.add_argument("-o", "--out", dest="outdir", type=str,
                  required=False, default="./", help="directory of output files")

# list output  name
#parser.add_argument("-l", "--list", dest="outlist", type=str,
#                  required=False, help="write list to OUT_LIST", metavar="OUT_LIST")

# base difference in snp position
parser.add_argument("-d", "--diff", dest="baseDiff", type=int,
                  required=False, default=10000, help="int of number of bases between snps")


parser.add_argument("-n", "--num", dest="snpCount", type=int,
                  required=False, default=500, help="number of snps to return")

# define args
args = parser.parse_args()



###########################################
############## Function(s) ################
###########################################

def getPositions(bamIn, chr):
	""" get snp positions from chromosome input """

	bamfile = pysam.AlignmentFile(bamIn, "rb")

	store = []
	print(chr)
	for read in bamfile.fetch(chr):
		store.extend(read.get_reference_positions())
	return np.unique(store)

def fuzzyDistance(diff, in_list):

	""" fuzzy match number to the next HIGHEST value """

	counter = 0
	store = []
	for i in in_list:
		if i > counter + diff or counter == 0:
			counter = i
			store.append(i) 
	return store


### save a text file without a new line at the end
def saveTxt(dirfile, listToSave, sep='\n'):
	with open(dirfile, 'w') as f:
		for num, val in enumerate(listToSave):
			if num == len(listToSave) - 1:
				f.write(val)
			else:
				f.write(val + sep)


def snpList(bamIn, pickleOut, baseDiff):

	""" return snp list from in bam file with 
		baseDiff number of bases apart"""

	chrPos = ["chr"+ str(i) for i in range(1,32)]
	chrPos.extend(["chrX", "chrM"]) # 32 chrX and 33 chrM

	positionList = [getPositions(bamIn, i) for i in chrPos]

	# fuzzy match next closest position
	fuzzyList = [fuzzyDistance(baseDiff, i) for i in positionList]

	# convert to strings
	fuzzyName = []
	for num, elem in enumerate(chrPos):
		print(str(num) + "/" + str(len(chrPos)))
		tmp = [elem + ":" + str(i) for i in fuzzyList[num]]
		fuzzyName.extend(tmp)

	# pickle 
	f = open(pickleOut,'wb') ## note the b: accept binary files
	pickle.dump(fuzzyName, f)
	f.close()



def subSample(pickleIn, listOut, snpCount):
	
	""" return a subsampled list of length: snpCount """
	
	f = open(pickleIn, 'rb')
	snps = pickle.load(f)
	f.close()
	
	np.random.seed(12345)
	subsample = np.random.choice(snps, snpCount, replace=False)
	subsample = natsorted(subsample)
	
	saveTxt(listOut, subsample)


def chrmSep(snpList, outDir):
	""" split snp lists by chromosome """
	f = open(snpList, "r")
	snps = f.read().splitlines()

	snps_split = [i.split(":") for i in snps]
	df = pd.DataFrame(snps_split)
	grouped = df.groupby(0)[1].apply(list)
	df_grouped = pd.DataFrame(grouped)


	# save files #
	for index, row in df_grouped.iterrows():
		outname = outDir + "." + row.name + ".list"
		content = [row.name + ":" + i for i in row[1]]
		saveTxt(outname, content)
	
	return 0


###########################################
############### Wraggling #################
###########################################
# python ancestry/selectSites.py -c snps -i data/processed_sequences/new.rg.bam -o data/ancestry/snp.all -d 10000

# python ancestry/selectSites.py -c subsample -i data/ancestry/snp.all.p -o data/ancestry/snp -n 5000

# python ancestry/selectSites.py -c split -i data/ancestry/snp.list -o data/ancestry/snp.chr/snp

if args.command == "snps":
	### generate pickle file ###
	snpList(bamIn=args.infile, \
		pickleOut=args.outdir + ".p", baseDiff=args.baseDiff)

elif args.command == "subsample":

	#### subsample ###
	subSample(pickleIn=args.infile, \
		listOut=args.outdir+".list", snpCount=args.snpCount)

elif args.command == "split":
	### split chromosomes ###
	chrmSep(snpList=args.infile, outDir=args.outdir)

else:
	print("not found")







