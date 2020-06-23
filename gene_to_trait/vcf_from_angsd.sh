#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-12
# Last Modified: 2020-06-18
# Desc: 
module load bcftools
module load anaconda3/personal

bcftools view gl.out.bcf > gl.out.vcf
#bcftools view $BCF | bcftools sort - > gl.out.vcf
cat gl.out.vcf | cut -f1-9,54 | gzip > target.vcf.gz

OUT=$DATA_DIR/annot.target

gunzip target.vcf.gz

# works online, but not on hpc
vep --biotype --check_existing --sift b --species equus_caballus --transcript_version --cache --tab --plugin Phenotypes,dir=$HOME/.vep/Plugins/,phenotype_feature=1 --input_file target.vcf --output_file $OUT