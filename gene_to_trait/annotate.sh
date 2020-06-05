#! /bin/bash

module load anaconda3/personal

DATA_DIR=data/gene_variants/annotate/
REF=data/processed_sequences/ref_genome/EquCab3.fna
BAM_FILE=data/processed_sequences/benson/final.bam

VCF=$DATA_DIR/final.chr3.raw.vcf
VCF_AN=$DATA_DIR/final.chr3.raw.annot.vcf
VCF_AN_FILTER=$DATA_DIR/final.chr3.filter.annot.vcf

#SNP_EFF=dependancies/snpEff/snpEff.jar

# snp calling of individual sample

#chr3:60000000-80200000
#chr3:77731743-77735488
#chr3:77731730-77735500
#chr3:79504300-79593715 # variable regions causing traits 

#samtools mpileup -uf $REF $BAM_FILE -r chr3:79504300-79593715 | bcftools call -m > $VCF
# bcftools -v variants only, remove for all sites



	# perhaps some QC here


	# annotate
	#java -Xmx4g -jar $SNP_EFF EquCab2.86 $VCF -v > $VCF_AN

	# move misc files
	#mv snpEff_genes.txt snpEff_summary.html $DATA_DIR

	#  filter variants not described
	#java -Xmx4g -jar dependancies/snpEff/SnpSift.jar filter -f $VCF_AN "! exists ID" > $VCF_AN_FILTER


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

	#./vep --appris --biotype --buffer_size 5000 \
	#	--check_existing --distance 5000 --mane \
	#	--plugin Phenotypes,dir=[path_to]/,phenotype_feature=1,exclude_sources=COSMIC&HGMD-PUBLIC&Cancer_Gene_Census \
	#	--sift b --species equus_caballus --symbol --transcript_version \
	#	--tsl --cache \
	#	--input_file [input_data] \
	#	--output_file [output_file]


### new command ###
#--plugin Phenotypes,dir=[path_to]/,phenotype_feature=1,exclude_sources=COSMIC&HGMD-PUBLIC&Cancer_Gene_Census \

# use json output for full phenotype effect, tab is fine for testing
	# --json 
	# --tab
cd $EPHEMERAL/sandbox

# working
vep --biotype --check_existing --sift b --species equus_caballus --transcript_version --cache --tab --plugin Phenotypes,dir=$HOME/.vep/Plugins/,phenotype_feature=1 --input_file final.chr3.raw.vcf --output_file annot.out

vep --appris --biotype --buffer_size 5000 --check_existing --distance 5000 --mane --plugin Phenotypes,dir=$HOME/.vep/Plugins/,phenotype_feature=1 --sift b --species equus_caballus --symbol --transcript_version --tsl --cache --input_file out.bref3.vcf.gz --output_file annot.out
