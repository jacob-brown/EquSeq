#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-22
# Last Modified: 2020-07-10



""" extract sites from geno file 
	and compare to info list """

###########################################
################# Modules #################
###########################################

import csv
import pandas as pd
import numpy as np
import gzip

###########################################
############## Function(s) ################
###########################################

###########################################
######### Input(s) and Parameters #########
###########################################


###########################################
######### Input(s) and Parameters #########
###########################################

# probabilities
#snp_info = pd.read_csv("data/gene_variants/snp.list.info.csv")
snp_info = pd.read_csv("data/gene_variants/snp.list.alt.info.csv")
fi = gzip.open('results/gene_to_trait/gl.out.geno.gz', 'rb')
geno_str = fi.readlines()
geno = [i.decode().replace("\n", "").split("\t") for i in geno_str]
devna = [i.remove("") for i in geno]

# get header from beagle file
fbeag = gzip.open('results/gene_to_trait/gl.out.beagle.gz', 'rb')
beagle_str = fbeag.readlines()
beagle_header = beagle_str[0].decode().replace("\n", "").split("\t")
del beagle_header[0:3]
beagle_header = ["chr", "pos"] + beagle_header

# convert to df and add header
geno_df = pd.DataFrame(geno)
geno_df.columns = beagle_header
# make marker string
marker_ar_geno = np.array(geno_df["chr"].astype(str) + "_" + geno_df["pos"].astype(str))
geno_df["marker"] = marker_ar_geno

# allele freq
fmafs = gzip.open('results/gene_to_trait/gl.out.mafs.gz', 'rb')
mafs_str = fmafs.readlines()
mafs = [i.decode().replace("\n", "").split("\t") for i in mafs_str]
mafs_df = pd.DataFrame(mafs[1:])
mafs_df.columns = mafs[0]

# get positions and count data
	# each line of the counts corresponds to the pos file
fpos = gzip.open('results/gene_to_trait/gl.out.pos.gz', 'rb')
poscount_str = fpos.readlines()
poscount = [i.decode().replace("\n", "").split("\t") for i in poscount_str]
poscount_marker = [i[0] + "_" + i[1] for i in poscount]
del poscount_marker[0]

fcount = gzip.open('results/gene_to_trait/gl.out.counts.gz', 'rb')
count_str = fcount.readlines()
count = [i.decode().replace("\n", "").split("\t") for i in count_str]
devna = [i.remove("") for i in count]
count_df = pd.DataFrame(count[1:])
count_header = [i.replace("TotDepth", "").capitalize() for i in count[0]]
count_df.columns = count_header
count_df.insert(loc=0, column='marker', value=np.array(poscount_marker))


# join allele freqs - markers and alleles
marker_ar_maf = np.array(mafs_df["chromo"].astype(str) + "_" + mafs_df["position"].astype(str))
mafs_df["marker"] = marker_ar_maf
mafs_df = mafs_df[["marker", "major", "minor"]]
mafs_df.columns = ["marker", "allele1", "allele2"]
data_df = pd.merge(mafs_df, geno_df, how="inner", on="marker")
data_df = data_df.drop(["chr", 'pos'], axis = 1)

# substitutions only
subs = snp_info[snp_info["mutation"] == "sub"]
subs = subs.copy() # not a proper copy yet
marker_ar = np.array("chr" + subs["chr"].astype(str) + "_" + subs["start_bp"].astype(str))
subs["marker"] = marker_ar
subs_trim = subs[["marker", "phen", "ref", "new"]]

# GP only subs
geno_subs = pd.merge(subs_trim, data_df, how="inner", on="marker")

# counts only subs
subs_count_df = count_df[count_df.marker.isin(np.array(geno_subs.marker))]

geno_subs.to_csv("results/gene_to_trait/genosub.csv")
subs_count_df.to_csv("results/gene_to_trait/markercount.csv")