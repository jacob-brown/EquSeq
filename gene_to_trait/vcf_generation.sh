#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-12
# Last Modified: 2020-06-18
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


# local testing
REF=data/processed_sequences/ref_genome/EquCab3.fna
BAM_FILE=data/processed_sequences/benson/final.bam


for i in {0..9}
do
	SNP_FILE=data/gene_variants/trait.snps/trait.snp.$i.bed
	VCF=data/gene_variants/vcf/trait.$i
	
	samtools mpileup -uf $REF $BAM_FILE -l $SNP_FILE | bcftools call -m | vcftools --vcf - --minQ 15 --out $VCF --recode --recode-INFO-all &
done; wait

#VCF_FILES=($(ls -d data/gene_variants/vcf/*.vcf))
ls -d data/gene_variants/vcf/*.vcf > data/gene_variants/vcf/vcf.list
bcftools concat -f data/gene_variants/vcf/vcf.list -O v > data/gene_variants/vcf/merged.vcf

gzip data/gene_variants/vcf/merged.vcf

scp data/gene_variants/vcf/merged.vcf.gz jb1919@login.cx1.hpc.ic.ac.uk:/rds/general/user/jb1919/ephemeral/gene_to_trait/