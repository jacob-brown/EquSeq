#! /bin/bash


# 1. generate GLs then run the following

#### VERSION 3 - Direct GL imputing ######

	# V3 beagle - GLs can be used in current format
	wget http://faculty.washington.edu/browning/beagle/beagle.jar
	# https://faculty.washington.edu/browning/beagle/b3.html

	#java -jar beagle.18May20.d20.jar gt=snps.chr13.raw.vcf.gz out=new_out chrom=chr13:761-761000

	#java -Xmx16g -jar beagle.27Jan18.7e1.jar gt=sortieMSL.recode.vcf out=result ref=chr22.1kg.vcf.gz chrom=22 impute=true

	#zcat < ALL.merged.beagle.gz | head -n500 > tmp.beagle
	#gzip tmp.beagle

	zcat < chr3.gl.out.beagle.gz | head > tmp.beagle
	gzip tmp.beagle

	java -Xmx6g -jar beagle.jar like=tmp.beagle.gz out=imp

	# remember to impute priot to using the sliding window and subsetting!
	zcat < imp.tmp.beagle.gz.phased.gz
	zcat < imp.tmp.beagle.gz.dose.gz
	zcat < imp.tmp.beagle.gz.gprobs.gz
	cat imp.tmp.beagle.gz.r2


###### ALT - new beagle, no GLs #######
cd test_beagle
wget http://faculty.washington.edu/browning/beagle/beagle.18May20.d20.jar

wget http://faculty.washington.edu/browning/beagle/bref3.18May20.d20.jar

wget http://faculty.washington.edu/browning/beagle/test.18May20.d20.vcf.gz


### VCF file ###

DATA_DIR=data/gene_variants/annotate/
REF=data/processed_sequences/ref_genome/EquCab3.fna
BAM_FILE=data/processed_sequences/benson/final.bam

VCF=sandbox/raw.vcf

# snp calling of individual sample

#chr3:79504300-79593715 # variable regions causing traits
# -b $EPHEMERAL/ancestry/bam.list

### on hpc ###
module load samtools/1.3.1
module load bcftools/1.3.1

samtools mpileup -uf $EPHEMERAL/ref_genome/EquCab3.fna -b $EPHEMERAL/ancestry/bam.list -r chr3:79504300-79593715 | bcftools call -m > raw.vcf
gzip raw.vcf

# scp to local

# test vcf is compatible
java -jar beagle.18May20.d20.jar gt=raw.vcf.gz out=imp.vcf

rm imp*
# vcf dimensions
zcat < raw.vcf.gz | tail -1 | tr '\t' '\n' | wc -l
# vcf is 54 columns wide
	# 54 is Benson
	# rest can be used as a ref

# split benson from the others
#zcat < raw.vcf.gz | cut -f1-53 > reference.vcf 
zcat < raw.vcf.gz | cut -f1-53 | tr '/' '|' | gzip > reference.raw.vcf.gz
zcat < reference.raw.vcf.gz | tail

# filter ungenotyped sites
bcftools filter --exclude "GT='.|.'" reference.raw.vcf.gz | gzip > reference.vcf.gz

#zcat < raw.vcf.gz | cut -f1-9 -f54 > target.vcf
zcat < raw.vcf.gz | cut -f1-9,54 | gzip > target.vcf.gz
zcat < target.vcf.gz | tail

# test with ref
java -jar beagle.18May20.d20.jar ref=reference.vcf.gz  gt=target.vcf.gz out=imp.vcf

# make bref
java -jar bref3.18May20.d20.jar reference.vcf.gz > reference.bref3

# run with bref
java -jar beagle.18May20.d20.jar ref=reference.bref3 gt=target.vcf.gz out=out.bref3


zcat < out.bref3.vcf.gz | tail -20
zcat < target.vcf.gz | tail -20


#### now VEP??? ####



