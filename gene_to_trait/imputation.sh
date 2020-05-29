#! /bin/bash


module load shapeit

 shapeit --input-bed binary_fileset.bed binary_fileset.bim binary_fileset.fam \
        --input-map genetic_map.txt \
        --output-max binary_fileset.phased.haps binary_fileset.phased.sample 


# 0.1 convert current file format
	# beagle to any format, then onto plink
	# when generating GLs use: 
	-doPlink 2 -doGeno -4
	#to generate tfam anf tpeg
# 0.1 alt use plink on vcf file


# remove multialleics
bcftools view --max-alleles 2 --exclude-types indels snps.chr13.raw.vcf.gz -o new.vcf.gz -Oz
# generate bed bim fam
./plink --vcf new.vcf.gz --maf 0.02 --make-bed --out binary_fileset 

# 0.2 genetic map
    # 1. The physical position (bp) [integer]
    # 2. The recombination rate (cM/Mb) [float]
    # 3. The genetic position (cM) [float]
    # pposition rrate gposition
	# 72765 0.12455 0.00000
	# 94172 0.12458 0.00266
	# 94426 0.12461 0.00269
	# 95949 0.12461 0.00288
	# 98087 0.12460 0.00315


# 0.2 make a reference panel

# 1. align gwas data

# 2. pre-phasing
shapeit -B binary_fileset -O gwas.phased 




# 3. impute
impute2 -use_prephased_g \
        -known_haps_g gwas.phased.haps \
        -h reference.haplotypes.gz \
        -l reference.legend.gz \
        -m genetic_map.txt \
        -int 9.1e6 9.6e6 \
        -Ne 20000 \
        -o gwas.from9.1e6_to9.6e6.imputed 


# ALT. beagle 
wget http://faculty.washington.edu/browning/beagle/beagle.jar


# alt beagle new
wget http://faculty.washington.edu/browning/beagle/beagle.18May20.d20.jar

java -jar beagle.18May20.d20.jar gt=snps.chr13.raw.vcf.gz out=new_out chrom=chr13:761-761000










