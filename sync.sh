#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-05-05
# Last Modified: 2020-07-06
# Desc: sync local files with those on HPC ignoring those in gitignore 

LOCAL_DIR=./
REMOTE_DIR=jb1919@login.cx1.hpc.ic.ac.uk:/rds/general/user/jb1919/home/genomics/EquSeq/


rsync -av --delete --include='data/ancestry/*' \
		--include='data/cleaned_data'\
		--include='data/*merge.csv'\
		--include='data/gene_variants/trait.snps/trait.snp.mend.list'\
		--include='checkEphemeral.py'\
		--include='checkEphemeral.sh'\
		--include='data/ancestry/bam_list_grps/*'\
		--include='data/snp_calling_list/*'\
		--exclude-from='.gitignore' \
		--include='data/scripts/copy.sh'\
		--exclude='.git/' \
		$LOCAL_DIR $REMOTE_DIR



