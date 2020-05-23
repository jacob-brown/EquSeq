#!/bin/bash
#PBS -lwalltime=24:00:00
#PBS -lselect=1:ncpus=1:mem=1gb

# Desc: Get files from ncbi sra

	# for sra_runs.txt
	# qsub -J 0-47 fastqDownload.sh

	# for sra_runs_to_do.txt
	# qsub -J 0-581 fastqDownload.sh
		# 0-10 test

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


# arg1 sra_list location
# arg2 fasta path out
python3 $CODE_DIR/setup/getFastq.py $CODE_DIR/data/cleaned_data/sra_runs.txt $RES_DIR

timer



echo '=================================='
echo -e "\nGet rest of Fasta files\n"


# arg1 sra_list location
# arg2 fasta path out
#python3 getFastq.py $EPHEMERAL/wgs_data/sra_runs_to_do.txt $EPHEMERAL/wgs_data/files_to_align/

timer

