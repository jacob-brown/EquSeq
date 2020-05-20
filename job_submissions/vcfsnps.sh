#! /bin/bash
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=4:mem=4gb

# qsub -J 0-33 vcfsnps.sh

module load anaconda3/personal

CODE_DIR=$HOME/genomics/EquSeq/
DIR=$EPHEMERAL/snp_calling/
#echo "${ALL_CHR[*]}"
#CHR="${ALL_CHR[1]}"
#CHR="${ALL_CHR[$PBS_ARRAY_INDEX]}"

echo "running..."

ALL_FILE=($(ls $DIR/*.vcf))
unset "ALL_FILE[10]"
#echo "${ALL_FILE[*]}"
#FILE="${ALL_FILE[$PBS_ARRAY_INDEX]}"


for FILE in "${ALL_FILE[@]}";do

	echo $FILE

	python3 $CODE_DIR/ancestry/selectSites.py -c snpsVCF \
		-i $FILE \
		-o $DIR/snp_lists/snp

done
