#!/bin/bash
#PBS -lwalltime=48:00:00
#PBS -lselect=1:ncpus=1:mem=1gb

# Desc: Get files from ncbi sra
	# wc -l 
	# qsub -J 0-629 fastqDownload.sh
		# qsub -J 0-3 fastqDownload.sh
		# qsub -J 3-629 fastqDownload.sh

CODE_DIR=$HOME/genomics/EquSeq/

# import unix functions
source $CODE_DIR/scripts/unix_functions.sh

DIR=$EPHEMERAL/wgs_data/
RES_DIR=$DIR/raw_files


echo '=================================='
echo -e "\nLoad modules\n"
module load anaconda3/personal



echo '=================================='
echo -e "\nGet Fasta files\n"

# 629 files - some require merging

# arg1 sra_list location
# arg2 fasta path out
python3 $CODE_DIR/setup/getFastq.py $CODE_DIR/data/cleaned_data/sra_runs_to_do.txt $RES_DIR
timer



#############################################
###### Reruns for failed downloads  ########
	# explore fastq_download_check.sh
	# check which codes have failed
	# 1. delete the partailly downloaded files
	# 2. rerun the above command 

# uncomment
#RAW_FILES=($(ls $RES_DIR))
#RERUNS=$CODE_DIR/data/cleaned_data/rerun.codes.txt

### remove the files - be careful! 
	# don't sub as a job, run from command line once

# uncomment
#for file in "${RAW_FILES[@]}"; do
#	BASE=($(echo $file | cut -d _ -f 1))
#	if grep -q $BASE $RERUNS; then
#		echo "removing: " $file
#		rm $file
#	fi
#done

### rerun on the issue list ###
	# comment out the proper script and run the following
#wc -l $RERUNS
	# qsub -J 0-66 fastqDownload.sh
#python3 $CODE_DIR/setup/getFastq.py $RERUNS $RES_DIR
#timer







#### MISC #####
# GB EU only
#python3 $CODE_DIR/setup/getFastq.py $CODE_DIR/data/cleaned_data/sra_runs.txt $RES_DIR

# 50 individuals


