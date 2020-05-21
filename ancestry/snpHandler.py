#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-05-04
# Last Modified: 2020-05-21



""" Generate sites for use in 
	genotype liklihood and PCA  """

###########################################
################# Modules #################
###########################################

import numpy as np
import time
from natsort import natsorted
import argparse
#import pandas as pd 
import os
import re
import subprocess


###########################################
################# Options #################
###########################################

parser = argparse.ArgumentParser(description=\
		'Generate snp position list on a fuzzy position of a bam file.')


# command 
#choose = ["snpsVCF", "snps", "subsample", "split", "snpsWG"]

choose = ["snps", "window", "subBeagle"]

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

# input list file
parser.add_argument("-l", "--listin", dest="listIn", type=str, required=False, \
						help="list input of strings, eg. chr:position")

# base difference in snp position
parser.add_argument("-d", "--diff", dest="baseDiff", type=str,
                  required=False, default=5, \
                  help="int kb bases between snps")


# define args
args = parser.parse_args()



###########################################
############## Function(s) ################
###########################################


##################################
############## misc ##############

### timer ###
def timer(state="S"):
	global startTimer

	"""" timer(), timer("e") to start and stop timer """
	if state.lower()=='s':
		
		print("timer started")
		startTimer = time.time() # start the timer from import	
	
	elif state.lower()=='e':

		duration = time.time() - startTimer
		duration = round(duration)
		string = "\n..........................\n"\
				"   Time elapsed: {} sec"\
				"\n..........................\n"\
				.format(duration)

		print(string)
	
	else:
		stop("state not found, s or e only.")


### save a text file without a new line at the end
def saveTxt(dirfile, listToSave, sep='\n'):
	with open(dirfile, 'w') as f:
		for num, val in enumerate(listToSave):
			if num == len(listToSave) - 1:
				f.write(val)
			else:
				f.write(val + sep)


###################################
######## Functions to call ########


### snp list from vcf file ### 
	# run for eeach chromosome
def vcfSNPs(vcfIn, listOut):
	timer()
	print("\nwriting snp list from vcf file: {}\n".format(vcfIn))
	chrome = os.path.basename(vcfIn).split(".")[1]
	fileOut =  listOut + "." + chrome + ".raw.list"
	
	command = "sed '/##/d' {} | cut -f1,2 ".format(vcfIn)

	print("converting to angsd readable")

	p = subprocess.Popen([command], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	out, er = p.communicate()
	out_string = out.decode()
	positions = out_string.replace("\t", ":").split("\n")
	positions.remove("")
	del positions[0]
	posSort = natsorted(positions)
	saveTxt(fileOut, posSort)
	timer("e")



### sliding window ###
def slidingWindow(beagleIn, listOut, distKB):

	distance = int(float(distKB) * 1000)

	### sites in beagle file ###
	no_ext = os.path.splitext(beagleIn)[0]
	command = "gunzip {} && sed '1d' {} | cut -f 1 | sed 's/_/:/' && gzip {}"\
						.format(beagleIn, no_ext, no_ext)

	p = subprocess.Popen([command], stdout=subprocess.PIPE, \
							stderr=subprocess.PIPE, shell=True)
	out, er = p.communicate()
	out_decode = out.decode()
	pos_str = out_decode.split("\n")
	pos_str.remove("")

	bePos = [i.split(":") for i in pos_str]

	chromo = np.unique([i[0] for i in bePos])

	if len(chromo) > 1:
		raise ValueError("chromosome number should now exceed 1. beagle file contains more.")

	pos_int = np.array([int(i[1]) for i in bePos])
	
	# bin data
	maxPos = max(pos_int)
	bins = np.array(range(1, maxPos, distance))
	bin_ind = np.digitize(pos_int, bins)

	# sample
	binned_pos_int = [pos_int[bin_ind == i] for i in range(1, len(bins))]

	sample = [str(np.random.choice(i, 1)[0]) for i in binned_pos_int if len(i) != 0]
	sample_str = [chromo[0] + "_" + i for i in sample]
	
	# write
	f = open(listOut, "w")
	f.write("\n".join(sample_str))
	f.close()




### subset beagle file ###
def subBeagle(beagleIn, beagleOut, posList):

	no_ext = os.path.splitext(beagleIn)[0]
	outFile = beagleOut + ".beagle"

	command_main = "gunzip {} && grep -Fw -f {} {} | cat sed -n '1p' {} - > {}"\
				.format(beagleIn, posList, no_ext, no_ext, outFile)

	devna = subprocess.Popen([command_main],  \
				stderr=subprocess.PIPE, shell=True)



###########################################
################## Main ###################
###########################################


####
## VCF - SNPs

### snps from vcf file ###
# python3 ancestry/snpHandler.py -c snps -i data/processed_sequences/snps/snps.chr11.raw.vcf -o data/ancestry/snp.chr/snp

### sliding window ###
# python3 ancestry/snpHandler.py -c window -i sandbox/snp.chr28.list.beagle.gz -o sandbox/chr28.list -d 13.5

### subset beagle from sliding window ###
# python3 ancestry/snpHandler.py -c subBeagle -i sandbox/snp.chr28.list.beagle.gz -o sandbox/snp.chr28.fi -l sandbox/chr28.list





if args.command == "snps":
	
	vcfSNPs(vcfIn=args.infile, listOut=args.outdir)

elif args.command == "window":

	slidingWindow(beagleIn=args.infile, \
					listOut = args.outdir, distKB=args.baseDiff)

elif args.command == "subBeagle":

	subBeagle(beagleIn=args.infile, \
				beagleOut = args.outdir, posList=args.listIn)

else:
	print("command not found")


