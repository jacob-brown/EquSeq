#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-05-18
# Last Modified: 2020-05-19

########### RE WRITE IF TIME PERMITS ###########
	# works but poorly designed
	# vcreate one file as I bind them together again in R


""" Convert report to a flat table """

###########################################
################# Modules #################
###########################################

#import numpy as np
import os
from natsort import natsorted

###########################################
############## Function(s) ################
###########################################


### save a text file without a new line at the end
def saveTaxa(filename, listToSave):
	with open(filename, 'w') as f:
	     for val in listToSave:
	     	f.write("\t".join(val) + "\n")


### read taxa ###
def readTaxRep(fileIn):
	""" read nested taxanomic report to a 
		nested list """
	data = []
	with open(file, "r") as f:
		reader = f.readlines()
		for i in reader:
			tmp = i.replace("\n", "").split("\t")
			data.append([c.strip() for c in tmp])
	return data

### filter rank ###
def rankRet(taxa, rank):
	
	""" return taxa for single rank level """
	
	#taxaIn = readTaxRep(fileIn)
	
	#head = ["perc", "nFrag", "nOnlyFrag", "rank", "ncbiTaxa", "name"]
	#unclassified = [i for i in taxaIn if i[3] == "U"]
	rank_taxa = [i for i in taxa if i[3] == rank]

	# update percentages and filter unwanted values
	#totReads = sum([float(i[1]) for i in rank_taxa])

	#for n, v in enumerate(rank_taxa):
	#	rank_taxa[n][0] = (float(v[1])/totReads)*100

	rank_taxa.sort(key = lambda x:float(x[0]), reverse=True)
	#fileName = os.path.basename(fileIn)
	#na = [i.append(fileName) for i in rank_taxa]
	#data = unclassified + rank_taxa 
	#data.insert(0, head)
	return rank_taxa

###########################################
######### Input(s) and Parameters #########
###########################################

report_dir = "results/oral_diversity/"
files = os.listdir(report_dir)
files_to_use = [report_dir + i for i in files if i.split(".")[1] == "kreport"]
files_to_use = natsorted(files_to_use)

###########################################
############### Wraggling #################
###########################################

domain = []
phylum = []
genus = []
species = []
unclassified = []
root = []
for file in files_to_use:

	taxa = readTaxRep(file)
	fileName = os.path.basename(file)

	tmp_d = rankRet(taxa, "D")
	na = [i.append(fileName) for i in tmp_d]
	domain = domain + tmp_d

	tmp_p = rankRet(taxa, "P")
	na = [i.append(fileName) for i in tmp_p]
	phylum = phylum + (tmp_p)
	
	tmp_g = rankRet(taxa, "G")
	na = [i.append(fileName) for i in tmp_g]
	genus = genus + (tmp_g)
	
	tmp_s = rankRet(taxa, "S")
	na = [i.append(fileName) for i in tmp_s]
	species = species + (tmp_s)

	tmp_u = rankRet(taxa, "U")
	na = [i.append(fileName) for i in tmp_u]
	unclassified = unclassified + (tmp_u)

	tmp_r = rankRet(taxa, "R")
	na = [i.append(fileName) for i in tmp_r]
	root = root + (tmp_r)
	
###########################################
############### Analysis ##################
###########################################

saveTaxa(report_dir+"/domain.txt", domain)
saveTaxa(report_dir+"/phylum.txt", phylum)
saveTaxa(report_dir+"/genus.txt", genus)
saveTaxa(report_dir+"/species.txt", species)
saveTaxa(report_dir+"/unclassified.txt", unclassified)
saveTaxa(report_dir+"/root.txt", root)

###########################################
############### Plotting ##################
###########################################