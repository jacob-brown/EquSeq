#! /bin/bash

######## VEP ###########
# job submitted online
	# appris - adds a flag, possibly only for humans
	# biotype - Adds the biotype of the transcript or regulatory feature
	# buffer_size - mem alocation default = 5000
	# check_existing - Checks for the existence of known variants
	# distance - upstream or downstream_gene_variant consequences, default = 5000
	# mane - mane flag, possibly only for humans
	# plugin - additional feats.
	# sift - amino acid change -b: prediction term and score
	# species
	# symbol -Adds the gene symbo
	# transcript_version - Add version numbers to Ensembl transcript identifiers 
	# tsl - transcript support level, possibly only for humans
	# cache

# use json output for full phenotype effect, tab is fine for testing
	# --json 
	# --tab

JOB_N=$1
DATA_DIR=$EPHEMERAL/gene_to_trait/
VCF=$DATA_DIR/trait.$JOB_N.recode.vcf
OUT=$DATA_DIR/annot.$JOB_N.raw.txt
#OUT=$DATA_DIR/annot.$JOB_N.raw.json --json


# if vcf is empty terminate session
LEN=($(grep -A 1 CHROM $VCF | wc -l ))
if (($LEN==1)); then 
	echo "vcf found no snps, no annoatation."
	echo "exiting." 
	exit
fi

# raw
vep --biotype --check_existing --sift b --species equus_caballus --transcript_version --cache --tab --plugin Phenotypes,dir=$HOME/.vep/Plugins/,phenotype_feature=1 --input_file $VCF --output_file $OUT


# working
# imputed
#vep --biotype --check_existing --sift b --species equus_caballus --transcript_version --cache --tab --plugin Phenotypes,dir=$HOME/.vep/Plugins/,phenotype_feature=1 --input_file out.bref3.vcf.gz --output_file annot.out
