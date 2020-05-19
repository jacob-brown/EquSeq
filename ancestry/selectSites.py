#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-05-04
# Last Modified: 2020-05-18



""" Generate sites for use in 
	genotype liklihood and PCA  """

###########################################
################# Modules #################
###########################################

import pysam
import numpy as np
import time
from natsort import natsorted
import argparse
import pandas as pd 
import os
import re


###########################################
################# Options #################
###########################################

parser = argparse.ArgumentParser(description=\
		'Generate snp position list on a fuzzy position of a bam file.')


# command 
choose = ["snpsVCF", "snps", "subsample", "split", "snpsWG"]
parser.add_argument("-c", "--command", dest="command", type=str,
					required=True, choices=choose, 
					help="choice function. snps: generate txt list of all sites; subsample: "\
						"the sites; chrmSplit: split outfiles by chromosome.")

# input bam file
parser.add_argument("-i", "--input", dest="infile", type=str,
					required=False, help="bam IN_BAM to have snps sampled from", \
					metavar="IN_BAM")

#  output directory
parser.add_argument("-o", "--out", dest="outdir", type=str,
                  required=False, default="./", help="directory of output files")

# base difference in snp position
parser.add_argument("-d", "--diff", dest="baseDiff", type=int,
                  required=False, default=10000, \
                  help="int of number of bases between snps")

# snp count
parser.add_argument("-n", "--num", dest="snpCount", type=int,
                  required=False, default=500, help="number of snps to return")

# all sites
parser.add_argument("-a", "--allsites", dest="allSNPs", action="store_true", \
					required=False, help="store TRUE for using all sites")

# define args
args = parser.parse_args()



###########################################
############## Function(s) ################
###########################################

def vcfSNPs(vcfIn, listOut):
	print("\nwriting snp list from vcf file\n")
	command = "sed  '/##/d' {} | awk 'NR>1' | cut -f1,2 > {}".format(vcfIn, listOut)
	os.system(command)




#def getPositions(bamIn, chr):
#	""" get snp positions from chromosome input """
#
#	bamfile = pysam.AlignmentFile(bamIn, "rb")
#
#	store = []
#	print(chr)
#	for read in bamfile.fetch(chr):
#		store.extend(read.get_reference_positions())
#	return np.unique(store)


### calculate the next number roughly (fuzzy) ### 
def fuzzyDistance(diff, in_list):

	""" fuzzy match number to the next HIGHEST value """

	# random start point within 0-diff
	store_fuzz = [np.random.choice(in_list[0:diff])]
	counter = store_fuzz[0] + diff # the next value

	for i in in_list:
		if i >= counter:
			counter = i + diff # update
			store_fuzz.append(i)

	return store_fuzz


### save a text file without a new line at the end
def saveTxt(dirfile, listToSave, sep='\n'):
	with open(dirfile, 'w') as f:
		for num, val in enumerate(listToSave):
			if num == len(listToSave) - 1:
				f.write(val)
			else:
				f.write(val + sep)


### all potential sites ###
def sitesWhole(diff, outDir):
	# base diff of report

	with open("data/EquCab2.0_assembly_report.txt", "r") as f:
		lines = [line.rstrip() for line in f]

	lengths = []
	chrom = []
	for elem in lines:
		if "#" not in elem and "scaffold" not in elem:
			tmp = elem.split()
			if tmp[7] == 'Primary':	
				lengths.append(int(tmp[9]))
				chrom.append(tmp[10])
			else:
				# mitochondria entry is shorter
				lengths.append(int(tmp[8]))
				chrom.append(tmp[9])

	

	### generate distances list ###
	length_diff = []
	for i in lengths:
		pos_Start = np.random.randint(1, diff, 1)[0]
		length_diff.append(np.array(range(pos_Start, lengths[0], diff)))
	outDir = "data/ancestry/snp.chr/snp"
	
	### write ###
	for n, c in enumerate(chrom):
		tmp_str = np.array([c + ":" + str(i) for i in length_diff[n]])
		filenam = outDir + "." + c + ".list"
		saveTxt(filenam, tmp_str)



	


### snp list of snps with a base difference apart ###
def snpList(snpList, txtOut, baseDiff):
	print("reading file")
	f = open(snpList,"r")
	
	pos = f.readlines()
	pos_update = [i.replace("\n","").split("\t") for i in pos]

	# unique chromosomes - remove scaffolds
	chromo = [i[0] for i in pos_update if not re.search('scaffold', i[0])]
	chromo_uniq = natsorted(np.unique(chromo))
	misc_remove = ["chrM", "chrX"]
	chromo_uniq = [i for i in chromo_uniq if i not in misc_remove]
	
	print("gathering...")
	store_snp = []
	for loc in chromo_uniq:
		print(loc)
		loc_pos = np.array([int(i[1]) for i in pos_update if i[0] == loc])
		loc_pos_uniq = np.sort(np.unique(loc_pos))
		store_snp.append([loc, loc_pos])


	# fuzzy match next closest position and convert to strings
	print("calculating distances with fuzzy increase.")
	fuzzyName = []
	for elem in store_snp:
		fuzzyList = fuzzyDistance(baseDiff, elem[1])
		tmp = [elem[0] + ":" + str(i) for i in fuzzyList]
		fuzzyName.extend(tmp)

	print("saving output.")
	
	saveTxt(txtOut, fuzzyName)
	return 0 
	
#def snpList(bamIn, txtOut, baseDiff):
#
#	""" return snp list from in bam file with 
#		baseDiff number of bases apart"""
#
#	chrPos = ["chr"+ str(i) for i in range(1,32)]
#	chrPos.extend(["chrX", "chrM"]) # 32 chrX and 33 chrM
#
#	print("getting position list from: ")
#
#	positionList = [getPositions(bamIn, i) for i in chrPos]
#
#
#	# fuzzy match next closest position
#	print("calculating distances with fuzzy increase.")
#	fuzzyList = [fuzzyDistance(baseDiff, i) for i in positionList]
#
#	# convert to strings
#	fuzzyName = []
#	for num, elem in enumerate(chrPos):
#		print(str(num) + "/" + str(len(chrPos)))
#		tmp = [elem + ":" + str(i) for i in fuzzyList[num]]
#		fuzzyName.extend(tmp)
#
#	# save snps txt 
#	saveTxt(txtOut, fuzzyName)


### subsample the snps ###
def subSample(txtIn, listOut, snpCount):
	
	""" return a subsampled list of length: snpCount """
	
	f = open(txtIn, "r")
	snps = f.read().splitlines()
	f.close()
	
	np.random.seed(12345)
	subsample = np.random.choice(snps, snpCount, replace=False)
	subsample = natsorted(subsample)
	
	saveTxt(listOut, subsample)


### generate seperate lists for each chromosome ###
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
# clear the wd
# rm data/ancestry/snp.all
# rm data/ancestry/snp.list
# rm data/ancestry/snp.chr/*


####
# Use genome report 
# 5kb apart
# python3 ancestry/selectSites.py -c snpsWG -o data/ancestry/snp.chr/snp -d 5000



####
## VCF - SNPs

### snps from vcf file ###
# python3 ancestry/selectSites.py -c snpsVCF -i data/processed_sequences/raw_variants.vcf -o data/ancestry/snp.vcf

# 5kb apart
# python3 ancestry/selectSites.py -c snps -i data/ancestry/snp.vcf.raw.list -o data/ancestry/snp.apart.list -d 5000


### split ###
# python3 ancestry/selectSites.py -c split -i data/ancestry/snp.apart.list -o data/ancestry/snp.chr/snp


### subset ###
# all
# python3 ancestry/selectSites.py -c subsample -i data/ancestry/snp.apart.list -o data/ancestry/snp -a


#### main ####
if args.command == "snps":
	### generate txt file ###
	snpList(snpList=args.infile, \
		txtOut=args.outdir, baseDiff=args.baseDiff)
elif args.command == "snpsWG":
	sitesWhole(diff=args.baseDiff, outDir=args.outdir)

elif args.command == "snpsVCF":
	vcfSNPs(vcfIn=args.infile, listOut=args.outdir+".raw.list")

elif args.command == "subsample":

	#### subsample ###
	# no of snps to use 
	if args.allSNPs:
		# all sites
		f = open(args.infile, "r")
		snps = len(f.readlines())
	else:
		snps=args.snpCount
	
	subSample(txtIn=args.infile, \
		listOut=args.outdir+".list", snpCount=snps)

elif args.command == "split":
	### split chromosomes ###
	chrmSep(snpList=args.infile, outDir=args.outdir)

else:
	print("not found")







