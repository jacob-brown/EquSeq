#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-22
# Last Modified: 2020-07-03



""" extract sites from beagle file 
	and compare to info list """

###########################################
################# Modules #################
###########################################

import csv
import pandas as pd
import numpy as np

###########################################
############## Function(s) ################
###########################################

###########################################
######### Input(s) and Parameters #########
###########################################

#snp_info = pd.read_csv("data/gene_variants/snp.list.info.csv")
snp_info = pd.read_csv("data/gene_variants/snp.list.alt.info.csv")
fi = open("results/gene_to_trait/gl.out.rn.beagle", "r")
beagle_str = fi.readlines()
beagle = [i.replace("\n", "").split("\t") for i in beagle_str]
header = beagle[0]
del beagle[0]
beagle_df = pd.DataFrame(beagle)
beagle_df.columns = header

###########################################
############### Wraggling #################
###########################################

# substitutions only
subs = snp_info[snp_info["mutation"] == "sub"]
subs = subs.copy() # not a proper copy yet
marker_ar = np.array("chr" + subs["chr"].astype(str) + "_" + subs["start_bp"].astype(str))
subs["marker"] = marker_ar
subs_trim = subs[["marker", "phen", "ref", "new"]]
beagle_subs = pd.merge(subs_trim, beagle_df, how="inner")

# remove issue sites  - unknown causes, investigate  later

beagle_subs.to_csv("results/gene_to_trait/beaglesub.csv")
