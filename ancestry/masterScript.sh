#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-01-21
# Last Modified: 2020-04-20
# Desc: Script that runs everything. Should be run if data is added/removed.

### clean all the old data ###
rm -r data/cleaned_data/*

### make new directory for info tables ###
mkdir data/cleaned_data/infotables_update

### run all the gathering scripts ###

# BioProject codes from supplementary materials is added
	# expanding and updating the master list
python3 ancestry/expandProject.py

# clean data and join supplementary tables
python3 ancestry/cleaningJoining.py

# generate a master list and summarise the run info data
python3 ancestry/summarise.py

# extract breeds of interest from list  - potentially tmp script
python3 ancestry/breeds.py

# which other files do we need?
python3 ancestry/getAllFastq.py