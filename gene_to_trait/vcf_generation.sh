#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-12
# Last Modified: 2020-06-12
# Desc: 

SNP_FILE=$1
JOB_N=$2

DATA_DIR=$EPHEMERAL/gene_to_trait/
REF=$EPHEMERAL/ref_genome/EquCab3.fna
BAM_FILE=$EPHEMERAL/novel_data/merged/final.bam

VCF=$DATA_DIR/trait.$JOB_N

# snp calling of individual sample

#chr3:60000000-80200000
#chr3:77731743-77735488
#chr3:77731730-77735500
#chr3:79504300-79593715 # variable regions causing traits 


samtools mpileup -uf $REF $BAM_FILE -l $SNP_FILE | bcftools call -m > $VCF.vcf

# filter bad sites
vcftools --vcf $VCF.vcf --minQ 15 --out $VCF --recode --recode-INFO-all

