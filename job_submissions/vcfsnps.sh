#! /bin/bash
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=4:mem=4gb

# qsub -J 0-33 vcfsnps.sh

module load anaconda3/personal

CODE_DIR=$HOME/genomics/EquSeq/
DIR=$EPHEMERAL/snp_calling/


echo "running..."

ALL_FILE=($(ls $DIR/*.vcf))
#echo "${ALL_FILE[*]}"
FILE="${ALL_FILE[$PBS_ARRAY_INDEX]}"


for FILE in "${ALL_FILE[@]}";do

	echo $FILE

	python3 $CODE_DIR/ancestry/snpHandler.py -c snpsVCF \
		-i $FILE \
		-o $DIR/snp_lists/snp

done

