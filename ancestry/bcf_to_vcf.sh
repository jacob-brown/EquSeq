#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-07-06
# Last Modified: 2020-08-20
# Desc: concat bcf files and convert to a single vcf file

module load bcftools/1.3.1; bcftools concat -O v -o merge.vcf *.bcf
