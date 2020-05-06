#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-05-05
# Last Modified: 2020-05-06
# Desc: sync local files with those on HPC ignoring those in gitignore 

LOCAL_DIR=./
REMOTE_DIR=jb1919@login.cx1.hpc.ic.ac.uk:/rds/general/user/jb1919/home/genomics/HorseGenomics/

rsync -av --delete --include='data/ancestry/*' \
		--exclude-from='.gitignore' \
		--exclude='.git/' \
		$LOCAL_DIR $REMOTE_DIR


