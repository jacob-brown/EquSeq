#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-05-04
# Last Modified: 2020-05-27



""" Functions for generate sites for use in 
	genotype liklihoods, and handling 
	beagle files """

###########################################
################# Modules #################
###########################################

import numpy as np
import time
from natsort import natsorted
import argparse
import os
import re
import subprocess
import time


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
parser.add_argument("-d", "--diff", dest="baseDiff", type=float,
                  required=False, default=5, \
                  help="float kb bases between snps")

# number of files of split snps by (bins)
parser.add_argument("-b", "--bins", dest="nbins", type=int,
                  required=False, default=10, \
                  help="number of bins to split snps by, saves a file for each bin.")

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
		
		print("\ntimer started...\n")
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
	# run for each chromosome
def vcfSNPs(vcfDir, listOut, nBins):
	timer()
	
	print("\nreading snps from vcf\n")

	command = "sed '/##/d' {}/*.vcf | cut -f1,2".format(vcfDir)

	p = subprocess.Popen([command], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	out, er = p.communicate()
	out_string = out.decode()
	positions = out_string.replace("\t", ":").split("\n")
	positions.remove("")
	header = positions[0]
	pos = np.array([i for i in positions if i != header])

	# split list into groups to be saved
	print("binning")
	split = np.array_split(pos, nBins)	
	
	for n, elem in enumerate(split):
		fileOut = listOut + "." + str(n) + ".list"
		saveTxt(fileOut, elem)

	print("done.")
	
	timer("e")

### snp window ###
def snpWindow(beagleIn, listOut, distKB):

	distance = int(float(distKB) * 1000)

	### sites in beagle file ###
	print("extracting all sites")
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
	chromo = [i[0] for i in bePos]

	pos_int = np.array([int(i[1]) for i in bePos])
	
	# bin data
	print("\nbinning sites\n")

	maxPos = max(pos_int)
	bins = np.array(range(1, maxPos, distance))
	bin_ind = np.digitize(pos_int, bins)

	# sample
	print("\nsampling\n")
	binned_pos_int = [pos_int[bin_ind == i] for i in range(1, len(bins))]

	sample = [str(np.random.choice(i, 1)[0]) for i in binned_pos_int if len(i) != 0]
	
	sample_str = [i[0] + "_" + i[1] for i in bePos if i[1] in sample]

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

	# tidy up and gzip
	gzip = lambda f : subprocess.Popen(["gzip", f])

	# allow for catchup, skipped when pauses are absent
	time.sleep(.1) 
	devna = gzip(outFile)
	time.sleep(.1)
	devna = gzip(no_ext)
	return 0


###########################################
################## Main ###################
###########################################


####
## VCF - SNPs

### snps from vcf file ###
# python3 ancestry/snpHandler.py -c snps -i data/processed_sequences/snps/ -o data/ancestry/snp.chr/snp -b 10

### snp window ###
# python3 ancestry/snpHandler.py -c window -i data/processed_sequences/beagle/snp.chr28.list.beagle.gz -o data/processed_sequences/beagle/chr28.list -d 5

### subset beagle from snp window ###
# python3 ancestry/snpHandler.py -c subBeagle -i data/processed_sequences/beagle/snp.chr28.list.beagle.gz -o data/processed_sequences/beagle/snp.chr28.fi -l data/processed_sequences/beagle/chr28.list




if args.command == "snps":
	
	vcfSNPs(vcfDir=args.infile, listOut=args.outdir, nBins=args.nbins)

elif args.command == "window":

	snpWindow(beagleIn=args.infile, \
					listOut = args.outdir, distKB=args.baseDiff)

elif args.command == "subBeagle":

	subBeagle(beagleIn=args.infile, \
				beagleOut = args.outdir, posList=args.listIn)

else:
	print("command not found")


