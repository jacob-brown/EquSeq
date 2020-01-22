#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-01-21
# Last Modified: 2020-01-21
# Desc: Script that runs everything. Should be run if data is added/removed.


# BioProject codes from supplementary materials is added
	# expanding and updating the master list
python3 expand_projects.py 

# Gather the infotables from NCBI
python3 infotable_gather.py 


