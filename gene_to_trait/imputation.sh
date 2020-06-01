#! /bin/bash


# currnet beagle new
wget http://faculty.washington.edu/browning/beagle/beagle.18May20.d20.jar

# ALT. beagle 
wget http://faculty.washington.edu/browning/beagle/beagle.jar
# https://faculty.washington.edu/browning/beagle/b3.html

java -jar beagle.18May20.d20.jar gt=snps.chr13.raw.vcf.gz out=new_out chrom=chr13:761-761000
java -Xmx16g -jar beagle.27Jan18.7e1.jar gt=sortieMSL.recode.vcf out=result ref=chr22.1kg.vcf.gz chrom=22 impute=true

zcat < ALL.merged.beagle.gz | head -n500 > tmp.beagle
gzip tmp.beagle

java -Xmx6g -jar beagle.jar like=tmp.beagle.gz out=res.out

# remember to impute priot to using the sliding window and subsetting!


zcat < res.out.tmp.beagle.gz.phased.gz | head -n2

zcat < res.out.tmp.beagle.gz.gprobs.gz






cat res.out.log