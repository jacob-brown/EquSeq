#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-02-10
# Last Modified: 2020-03-02


##############################################
######### fastqc report file control #########
##############################################
# script controls which files are run and the stored output


#iter <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX")) # iteration for cluster
iter <- 1 # for testing

# files - with paths
files <- list.files('/rds/general/user/jb1919/home/genomics/sequences/cdts-hk.genomics.cn/Clean/F19FTSEUHT1854-swab-horse-1A/',  full.names=TRUE) 

# path out
path_out <- '/rds/general/user/jb1919/home/genomics/results/report/'

# command
command <- sprintf('fastqc -d . -o %s %s', path_out, files[iter]) # path out and file name on position

# run command
system2(command)

