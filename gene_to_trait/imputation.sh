#! /bin/bash

#BEAGLE=$EPHEMERAL/dependencies/beagle.jar 
BEAGLE=../../dependancies/beagle.jar
#BREF=$EPHEMERAL/dependencies/bref.jar 
BREF=../../dependancies/bref.jar

#cd $EPHEMERAL/gene_to_trait
cd data/gene_variants/
BCF=gl.out.bcf

# 0. setup
#wget -O dependancies/beagle.jar https://faculty.washington.edu/browning/beagle/beagle.27Jan18.7e1.jar

#wget -O dependancies/bref.jar http://faculty.washington.edu/browning/beagle/bref.27Jan18.7e1.jar

# 1. generate GLs 
	# genotypeLiklihood.sh
# convert to vcf and zip
#module load bcftools

# to run on hpc 
source activate myenv

# 2. scp to local
#scp jb1919@login.cx1.hpc.ic.ac.uk:/rds/general/user/jb1919/ephemeral/gene_to_trait/gl.out.bcf data/gene_variants/



######################################
### ALT ###
	# BCF format from angsd
#module load bcftools/1.3.1
#bcftools view gl.out.bcf > gl.out.vcf
bcftools view $BCF | bcftools sort - > gl.out.vcf
gzip gl.out.vcf


# test vcf is compatible
# GT
#java -jar $BEAGLE gt=gl.out.vcf.gz out=imp.vcf

# GL
#java -jar $BEAGLE gl=gl.out.vcf.gz out=imp.vcf


### generat reference and target groups ###
# individual count
zcat < gl.out.vcf.gz | tail -1 | tr '\t' '\n' | wc -l

# reference
# 53 should change with the number of individuals
zcat < gl.out.vcf.gz | cut -f1-53 | tr '/' '|' | gzip > reference.vcf.gz
#zcat < reference.gl.out.vcf.gz | tail

# filter ungenotyped sites
#bcftools filter --exclude "GT='.|.'" reference.gl.out.vcf.gz | gzip > reference.vcf.gz

#zcat < reference.vcf.gz | tail -20

# target 
# 54 should change with the number of individuals
zcat < gl.out.vcf.gz | cut -f1-9,54 | gzip > target.vcf.gz
#zcat < target.vcf.gz | tail -20


# test with ref
#java -jar $BEAGLE ref=reference.vcf.gz  gl=target.vcf.gz out=imp.vcf

# make bref
java -jar $BREF reference.vcf.gz 

# run with bref
java -jar $BEAGLE ref=reference.bref gl=target.vcf.gz out=out.bref3

#zcat < out.bref3.vcf.gz | tail -20
conda deactivate


#4. now VEP
scp out.bref3.vcf.gz jb1919@login.cx1.hpc.ic.ac.uk:/rds/general/user/jb1919/ephemeral/gene_to_trait/

