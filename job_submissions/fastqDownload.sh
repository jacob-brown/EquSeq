#!/bin/bash
#PBS -lwalltime=24:00:00
#PBS -lselect=1:ncpus=1:mem=1gb

# Desc: Get files from ncbi sra

# for sra_runs.txt
# qsub -J 0-47 fastqDownload.sh


# for sra_runs_to_do.txt
# qsub -J 0-581 fastqDownload.sh
	# 0-10 test


# import unix functions
cp $HOME/genomics/code/unix_functions.sh $TMPDIR
source unix_functions.sh

DIR=$EPHEMERAL/all_data/
RES_DIR=$DIR/files/


echo '=================================='
echo -e "\nLoad modules\n"
module load anaconda3/personal



echo '=================================='
echo -e "\nGet Fasta files\n"

#cp $HOME/genomics/code/getFastq.py $TMPDIR

# arg1 sra_list location
# arg2 fasta path out
#python3 getFastq.py $EPHEMERAL/all_data/sra_runs.txt $EPHEMERAL/all_data/files/

#timer



echo '=================================='
echo -e "\nGet rest of Fasta files\n"

cp $HOME/genomics/code/getFastq.py $TMPDIR

# arg1 sra_list location
# arg2 fasta path out
python3 getFastq.py $EPHEMERAL/all_data/sra_runs_to_do.txt $EPHEMERAL/all_data/files_to_align/

timer

