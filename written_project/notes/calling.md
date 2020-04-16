# Calling 

* genotype calling
* allele frequency estimation
* variant (or SNP) calling

## packages
angsd 
https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-014-0356-4

https://github.com/mfumagalli/Copenhagen/blob/master/Files/day1.md
## General
ANGDS uses a base set of functions for filtering, then uses a variaty of other functions to calculate statistics. Eg. 

```
$NGS/angsd/angsd -b $DATA/EUR.bams -ref $REF -out Results/EUR \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 7 -setMaxDepth 30 -doCounts 1 \
	-GL 2 -doGlf 1

$NGS/angsd/angsd -glf Results/EUR.glf.gz -fai $REF.fai -nInd 10 -out Results/EUR \
	-doMajorMinor 1 -doGeno 3 -doPost 2 -doMaf 1

```

## Steps
generate .fai index file with samtools
ensure all correct headers are present - RG isn't

1) Data filtering and I/O
2) Genotype likelihoods

`/angsd -GL`
3) Genotype calling
4) SNP calling

