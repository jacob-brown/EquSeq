#! /bin/bash
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=8:mem=10gb

# qsub vcfsnps.sh

module load anaconda3/personal

CODE_DIR=$HOME/genomics/EquSeq/
DIR=$EPHEMERAL/snp_calling/


#### 
# generate the snp list from vcf files
	#-b X, qsub -J 0-X ancestry.sh

python3 $CODE_DIR/ancestry/snpHandler.py -c snps \
			-i $DIR -o $DIR/snp_lists/snp -b 400

# 100 bins would take ~96hr