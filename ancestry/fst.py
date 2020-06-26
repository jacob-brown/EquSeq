#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-26
# Last Modified: 2020-06-26



""" calculate pairwise Fst using realSFS """

###########################################
################# Modules #################
###########################################

import os
import numpy as np
import subprocess

###########################################
############## Function(s) ################
###########################################

def expandGridUniq(listIn):
	
	""" expand a list returning 
		only unique combinations """
	
	# create an "all combos" matrix
		# removing those that are matches
	combos_all = np.array([[x, y] for x in listIn for y in listIn if x != y])
	
	# remove combinations that have already occured
		# resort and then filter uniques only
	combos_sorted = [np.sort(i) for i in combos_all]
	combos_true = [list(y) for y in set([tuple(x) for x in combos_sorted])]
	return combos_true

# sub process wrapper
def subProWrap(command, returnList=True):
	""" wrapper for a subprocess command 
		returns list as default """
	p = subprocess.Popen([command], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	out, er = p.communicate()
	if(returnList):
		out_string = out.decode()
		files = out_string.replace("\t", ":").split("\n")
		files.remove("")
		return files


# between populations 2dSFS
def twoDSFS(pop1, pop2, realSFS= "/rds/general/user/jb1919/ephemeral/dependencies/angsd/misc/realSFS"):

	""" run realSFS 2dSFS """ 

	cmd_twoDSFS = "{realsfs} {pop1}.out.saf.idx "\
					"{pop2}.out.saf.idx -fold 1 "\
					"-P 4 2> /dev/null > {pop1}.{pop2}.sfs"\
					.format(realsfs=realSFS, pop1=pop1, pop2=pop2)

	print("saving: " + "{}.{}.sfs".format(pop1, pop2))

	subProWrap(cmd_twoDSFS, returnList=False)


	# this can be improved with multiple pops in one  !!!!!
def fstpbs(pop1, pop2, realSFS= "/rds/general/user/jb1919/ephemeral/dependencies/angsd/misc/realSFS"):
	
	""" generate Fst files """

	cmd_siteFst = "{realsfs} fst index {pop1}.out.saf.idx "\
					"{pop2}.out.saf.idx -sfs {pop1}.{pop2}.sfs "\
					"-fstout {pop1}.{pop2}.pbs -whichFST 1 &> /dev/null"\
					.format(realsfs=realSFS, pop1=pop1, pop2=pop2)

	print("saving: " + "{}.{}.pbs.fst.idx".format(pop1,pop2))

	subProWrap(cmd_siteFst, returnList=False)


def fstStat(pop1, pop2, realSFS= "/rds/general/user/jb1919/ephemeral/dependencies/angsd/misc/realSFS"):

	""" get pairwise fst. 
		returns [FST.Unweight, Fst.Weight] """

	cmd_pairFst = "{realsfs} fst stats {pop1}.{pop2}.pbs.fst.idx"\
					.format(realsfs=realSFS, pop1 = pop1, pop2=pop2)

	fst_str = subProWrap(cmd_pairFst, returnList=True)
	print(fst_str[0].split(":"))
	return fst_str[0].split(":")


###########################################
######### Input(s) and Parameters #########
###########################################

list_path = "/rds/general/user/jb1919/home/genomics/EquSeq/data/ancestry/bam_list_grps/"
files = os.listdir(list_path)
pops = [i.replace(".list", "") for i in files]

# parse this as arg
realSFS = "/rds/general/user/jb1919/ephemeral/dependencies/angsd/misc/realSFS"

###########################################
############### Wraggling #################
###########################################

combos = expandGridUniq(pops)


twoDSFS(pop1='deutschesreitpony', pop2='connemara')
fstpbs(pop1='deutschesreitpony', pop2='connemara')
fstStat(pop1='deutschesreitpony', pop2='connemara')




###########################################
############### Analysis ##################
###########################################



###########################################
############### Plotting ##################
###########################################