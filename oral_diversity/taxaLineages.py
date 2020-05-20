#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-05-18
# Last Modified: 2020-05-20



""" Generate ncbi taxa dataframe from Kraken2 reports and taxonkit  """

###########################################
################# Modules #################
###########################################

import os
from natsort import natsorted
import subprocess 
import numpy as np
import pandas as pd
import argparse

###########################################
################# Options #################
###########################################

parser = argparse.ArgumentParser(description=\
		'Generate ncbi taxa dataframe from Kraken2 reports and taxonkit.')

# directory to kraken2 reports
parser.add_argument("-d", "--dirIn", dest="dirin", type=str,
					required=True, help="directory to kraken files")

#  output directory
parser.add_argument("-o", "--out", dest="outdir", type=str,
                  required=False, default="./", help="directory of output files")

#  output directory
parser.add_argument("-f", "--format", dest="format", type=str,
                  required=False, default="df", \
                  help="""df: output tab-deliminated dataframe; 
                  			l: one line list of lineage""")

# define args
args = parser.parse_args()

###########################################
############## Function(s) ################
###########################################

### get lineage of taxa ###
def taxLineage(dirIn, dirOut, outFormat):
	
	""" return lineage of taxa from kreports.
		provide directory to reports and
		ensure taxonkit is installed """

	### generate a unique list of all taxa ID ###
	command_id = "awk FNR-1 {}/*.kreport | cut -f 5".format(dirIn)
	p_id = subprocess.Popen([command_id], stdout=subprocess.PIPE,\
							stderr=subprocess.PIPE, shell=True)
	out, er = p_id.communicate()
	out_str = out.decode()
	data = out_str.split("\n")
	data.remove("")
	taxa_ID = np.unique(data)
	taxa_save = "\n".join(taxa_ID)
	
	outFile = dirOut + "/taxa.id"
	
	print("\nsaving taxaID: {}\n".format(outFile))

	f = open(outFile, "w")
	f.write(taxa_save)
	f.close()

	### run taxonkit ###

	### data frame format
	if outFormat.lower() =="df":
		command_taxa = """taxonkit lineage {} | 
							taxonkit reformat -p '' -r '' | 
							csvtk -H -t cut -f 1,3 | 
							csvtk -H -t sep -f 2 -s ';' -R | 
							csvtk add-header -t -n 
								taxid,kindom,phylum,class,
								order,family,genus,species""".format(outFile)
	
	elif outFormat.lower() =="l":
		### listed format ###	
		command_taxa = "taxonkit lineage -r --show-status-code {}".format(outFile)
	else:
		stop("format unknown")


	p_taxa = subprocess.Popen([command_taxa], stdout=subprocess.PIPE,\
									stderr=subprocess.PIPE, shell=True)
	out_t, er_t = p_taxa.communicate()
	out_t_str = out_t.decode()
	data_lines = out_t_str.split("\n")
	data_lines.remove("")
	data = [i.split("\t") for i in data_lines]

	outDF = dirOut + "/taxa.df"

	print("\nsaving taxaID: {}\n".format(outFile))
	with open(outDF, 'w') as f:
		for val in data:
			f.write("\t".join(val) + "\n")

	return 0


###########################################
################# main ####################
###########################################

# python3 oral_diversity/taxaLineages.py -d results/oral_diversity/kraken_reports/ -o results/oral_diversity

# lineage data
taxLineage(dirIn=args.dirin, dirOut= args.outdir, outFormat=args.format)




