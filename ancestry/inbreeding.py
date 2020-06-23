#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-17
# Last Modified: 2020-06-23



"""  """

###########################################
################# Modules #################
###########################################

import numpy as np

###########################################
############## Function(s) ################
###########################################

inbreed = np.load("results/ancestry/ALL.PCA.inbreed.npy")
#admix = np.load("results/ancestry/ALL.PCA.admix.Q.npy")
f = open("results/ancestry/clusters", "r")
data_str = f.readlines()
clst = [i.replace("\n", "").split("\t")[2] for i in data_str]

data = []
for i in range(0,len(clst)):
	data.append([clst[i], inbreed[i]])


max(inbreed)
min(inbreed)

###########################################
######### Input(s) and Parameters #########
###########################################



###########################################
############### Wraggling #################
###########################################



###########################################
############### Analysis ##################
###########################################



###########################################
############### Plotting ##################
###########################################