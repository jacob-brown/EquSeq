#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-01-21
# Last Modified: 2020-01-28
# Desc: Script that runs everything. Should be run if data is added/removed.

### clean all the old data ###
rm ../data/infotables_update/* # updated runinfo tables

$ cat to_remove
/tmp/file1
/tmp/file2
/tmp/file3
$ rm $( cat to_remove )

rm ../data/papers_projects_update.csv
rm ../data/prj_ancient.txt
rm ../data/prj_modern.txt

### run all the gathering scripts ###

# BioProject codes from supplementary materials is added
	# expanding and updating the master list
python3 expand_projects.py 

# clean data and join supplementary tables
python3 cleaning_joining.py

# generate a master list and summarise the run info data
python3 summary.py


