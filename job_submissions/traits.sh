#! /bin/bash
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=5:mem=5gb

#qsub -J 0-9 traits.sh

module load samtools/1.3.1
module load bcftools
module load vcftools/0.1.14
module load anaconda3/personal

CODE_DIR=$HOME/genomics/EquSeq/
source $CODE_DIR/scripts/unix_functions.sh

#ALL_FILES=($( ls -v $CODE_DIR/data/gene_variants/trait.snps/*.bed ))
ALL_FILES=($( ls -v $CODE_DIR/data/gene_variants/trait.snps/*.tab ))
FILE=${ALL_FILES[$PBS_ARRAY_INDEX]} 


echo '=================================='
echo -e "\nGenerating VCF\n"
	# arg1: bed file
	# arg2: job number (naming)
sh $CODE_DIR/gene_to_trait/vcf_generation.sh $FILE $PBS_ARRAY_INDEX

timer


echo '=================================='
echo -e "\nAnnotating\n"

sh $CODE_DIR/gene_to_trait/annotate.sh $PBS_ARRAY_INDEX

timer

