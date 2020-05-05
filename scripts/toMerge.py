#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-04-16
# Last Modified: 2020-05-05



""" determine which fastq files need to be merged """



# python3 scripts/toMerge.py -b data/bam.list.all.txt -o data/to_merge.csv -i data/cleaned_data/info_all.csv -r Run -g BioSample

###########################################
################# Modules #################
###########################################

import csv
import os
import argparse
import numpy as np

###########################################
################# Options #################
###########################################

parser = argparse.ArgumentParser(description='Which files need merging.')

# input list
parser.add_argument("-b", "--bam.list", dest="bamlist", type=str,
					required=True, metavar="IN_LIST",
					help="IN_LIST.csv of unique RUN ID and grouping, eg. individual")

# output list
parser.add_argument("-i", "--in.file", dest="infile", type=str,
                 	 help="IN_FILE with information on RUN_ID and GROUP_ID",
                 	 required=True, metavar="IN_FILE")

# output list
parser.add_argument("-o", "--out.list", dest="outfile", type=str,
                 	 help="write OUT_LIST of RUN IDs to merge and the grouping",
                 	 required=True, metavar="OUT_LIST")

parser.add_argument("-r", "--runID", dest="runid", type=str,
                  required=True, help="RUN_ID name, should be unique.",
                  metavar="RUN_ID")

parser.add_argument("-g", "--groupby", dest="grp", type=str,
                  required=True, help="group by GROUP_ID variable name.",metavar="GROUP_ID")


# define args
args = parser.parse_args()

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

#bamlist = args.bamlist
#args.runid  # Run
#args.grp # BioSample
#args.outfile # '../data/cleaned_data/to_merge.csv'
#args.infile #'../data/cleaned_data/info_all.csv'

# test locally 
#echo ERR2179543.bam\n  SRR1769892.bam\n  SRR1769922.bam\n  SRR515202.bam\n  SRR515204.bam\n  SRR515206.bam\n  SRR515209.bam\n  SRR515212.bam\n  SRR515214.bam\n  SRR515216.bam\n ERR2731056.bam\n  SRR1769893.bam\n  SRR1790681.bam\n  SRR515203.bam\n  SRR515205.bam\n  SRR515208.bam\n  SRR515211.bam\n  SRR515213.bam\n  SRR515215.bam\n  ERR979130.bam\n ERR979131.bam\n ERR979133.bam\n ERR979134.bam\n  ERR979218.bam\n ERR979219.bam\n ERR979220.bam\n ERR979221.bam\n ERR979223.bam\n ERR979224.bam\n ERR979225.bam > ../sandbox/bam.list.txt

# python3 toMerge.py -b ../sandbox/bam.list.txt -o ../data/cleaned_data/to_merge.csv -i ../data/cleaned_data/info_all.csv -r Run -g BioSample

# HPC
# ls 


info_all = open_csv(args.infile)

bam_txt = open(args.bamlist, "r")
bam_txt = open("data/bam.list.all.txt", "r")
bams = [i.split()[0] for i in bam_txt.readlines()]

###########################################
############### Wraggling #################
###########################################


bamlist_no_ext = [i.split(os.extsep, 1)[0] for i in bams]

i_run = info_all[0].index(args.runid)
i_grp = info_all[0].index(args.grp)

# only files present in bamlist
present_only = [i for i in info_all if i[i_run] in bamlist_no_ext]

# unique codes - i.e. don't group or merge
groups = [i[i_grp] for i in present_only]

# return non unique grouping codes
not_unique_vals = [i for i in {*groups} if groups.count(i) > 1]

# and unique ones
unique_vals = [i for i in {*groups} if groups.count(i) == 1]

store = []
for val in not_unique_vals:
	tmp = [i[i_run] for i in present_only if val == i[i_grp]]
	tmp.insert(0, val)
	store.append(tmp)

# append a header
store.insert(0, [args.grp, args.runid + '*'])

# write to file
write_csv(store, args.outfile)

