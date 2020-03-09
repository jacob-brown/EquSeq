#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-01-21
# Last Modified: 2020-03-09
# Desc: Script that runs everything. Should be run if data is added/removed.

### clean all the old data ###
rm -r ../data/cleaned_data/*

### make new directory for info tables ###
mkdir ../data/cleaned_data/infotables_update

### run all the gathering scripts ###

# BioProject codes from supplementary materials is added
	# expanding and updating the master list
python3 getall_expand_projects.py 

# clean data and join supplementary tables
python3 getall_cleaning_joining.py

# generate a master list and summarise the run info data
python3 getall_summary.py


