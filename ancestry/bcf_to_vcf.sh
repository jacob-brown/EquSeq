#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-07-06
# Last Modified: 2020-07-07
# Desc: concat bcf files and convert to a single vcf file


module load bcftools/1.3.1; bcftools concat -O v -o merge.vcf *.bcf

#bcftools view snp.0.vcf.list.bcf

# scp jb1919@login.cx1.hpc.ic.ac.uk:/rds/general/user/jb1919/ephemeral/ancestry/all_bcf/merge.vcf.gz sandbox/snps/ALL.merge.vcf.gz