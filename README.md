# EquSeq

Submitted as part of the thesis:

**Title:** Direct to consumer genetics for non-human samples
**Author**: Jacob Brown
Master of Research at Imperial College London.
Submitted for the *MRes in Computational Methods in Ecology and Evolution*.


**N.B.** General usage and thesis replicability is sufficient, however, the tool is a work in progress, eg. many of the scripts are still hard-coded.

The scripts make use of 2 environments (both unix), HPC and local. Much of the analysis on raw sequences occured on the HPC space, whilst the local environmet was reserved for smaller data file handling and plotting. 

-----
# Files/Directories of particular interest
For thesis examinars and those undertaking future work, below are several important files/directories to take particular note of.

* `data/raw_data` - raw infotables from ncbi, compiled tables for scripts to pivot on. Only edit if more data is being added. 
* `data/metadata` - metadata tables from publications, usually breed-sample info 
* `data/cleaned_data` - scripts take rawdata, clean, assign breeds, and get dumped here
    * `info_all.csv` - Unique sequence runs/file codes with additional information
    * `info_individual_grouped.csv` - As above, but grouped at the individual sample level
    * `info_pop_grouped.csv` - As above, but grouped at the breed level
* `results/individs_used.csv` - assecion code's and corresponding breeds actually used

# Repo structure

```
.
├── ancestry
├── data
│   ├── cleaned_data
│   │   └── infotables_update
│   ├── gene_variants
│   │   ├── animalqtldb
│   │   ├── trait.snps
│   ├── metadata
│   ├── raw_data
│   │   ├── infotables_original
│   │   └── supplementary_data_from_studies
│   └── snp_calling_list
├── gene_to_trait
├── job_submissions
├── mapping
├── oral_diversity
├── results
│   ├── ancestry
│   │   ├── fstat
│   │   └── treemix
│   ├── gene_to_trait
│   └── oral_diversity
│       ├── kraken_reports
│       └── stats
├── scripts
└── setup

```
-----
# Dependencies 
Written in the unix environment.

## Languages
* Python (3.6.2)
* R (3.6.1)
* bash
* MySQL (scripts present but not utilised in the final report)

## Programs/Tools
* bcftools (1.3.1)
* Burrows-Wheeler Aligner (BWA) (0.7.8)
* fastqc (0.11.9)
* FASTX (0.0.14)
* KRAKEN (2.0)
* NGSadmix (v32)
* PCAngsd (0.985)
* picard (2.6.0)
* Plink (1.9)
* samtools (1.10)
* Taxonkit (0.6)
* Treemix (1.13)
* Vcftools (0.1.13)
* Plink2treemix.py

## Python Modules
* argparse
* copy
* csv
* gzip
* itertools
* math
* natsort
* numpy
* os
* pandas
* re
* requests
* scipy
* subprocess
* sys
* time

## R packages
* colorspace
* ggforce
* gghighlight
* ggplot2
* ggpubr
* gridExtra
* optparse
* pophelper
* RColorBrewer
* scales
* stringr
* tidyverse
* vegan

-----
# Usage/Notes

Always run the scripts from the parent directory `EquSeq`, unless specified.

`your_path` - notes the path you are using (ssh path for HPC if using) 

--------
# Setup
Programs requiring specific setup (or those that are slightly tricky) have specific scripts in the `setup` directory. Scripts here also are required for eg. pulling *fastq* files from ncbi into the desired location, and creating the `data/cleaned_data` contents. 



--------
# Data Preparation 
## Sequences
1. Novel sample
  - Mouth swab was taken from individual
  - low-pass WGS 
  - Illumina pair-end
  - 100bp long reads
  - BGI notes that qualities are in phred33 format
2. Published and publicly available
  - FASTQ files from ENA (European nucleotide archive)
  - WGS of horses, only those of known breeds
  - pair ended 
  - 38 breeds/mixed breed combinations (excluding novel)
  - 170 individuals (excluding novel)

## Publically available WGS data
1. NCBI BioProject code list and metadata on breeds, age of sample, etc. 
2. NCBI infotables downloaded for each NCBI BioProject (manually) as ftp leaves out some data
3. `masterScript.sh` is implemented joining supplementary materials, filtering, and summarising data 
  - results in lists runIDs to download from ENA/NCBI with the corresponding metadata
4. Data that passes has the raw sequence files downloaded from ENA

## Read Mapping - Novel Sample
`mapping/novelMapper.sh`
1. `fastqc` - check raw sequence quality
2. `fastx` - trim the last 10 bases of reads, due to poor kmer scores quality (our seq only)
3. `bwa mem` - align pair ended reads to reference genome
4. `samtools` - convert to bam format
5. `samtools` - sort read order
6. `samtools` - index
7. what were the flagstats? 
8. `samtools` - merge reads (our sample and if necessary)
9. `samtools` - index
10. `picard FixMateInformation` - check mate-pair information is in sync between each read and its mate pair
11. `picard MarkDuplicates` - remove duplicate reads
12. `samtools` - remove mapping quality below 20
13. our sample: `picard AddOrReplaceReadGroups` - add a read group
14. `samtools` - remove unmapped reads
15. `samtools` - re-index

## Read Mapping - WGS samples
0. `wgsMapper.sh` used to control the mapping, as pairing read identification is different.
    - number of jobs should equal the number of run codes (as some are paired and others not)
...

8. run `sh bamAncestry.sh -m` to generate the bam list on the hpc 
9. run `toMerge.py` to determine which samples need to be merged 
    * scp bam.list to local, investigate, and run (rsync one way compatible)
      - `scp your_path/ephemeral/ancestry/bam.list data/`
    * `python3 scripts/toMerge.py -b data/bam.list -o data/ -i data/cleaned_data/info_all.csv -r Run -g BioSample`
10. `rsync`
11. run `wgs_merge.sh` to merge the bam files
    * one job for each file (i.e. nrow of `to_merge.csv`)
    * also indexes bam files
12. move the merged files and files that don't require merging to the `final/` directory
    * hpc: `python ~/genomics/EquSeq/mapping/wgs_mover.py` - generate move list
    * hpc: `sh ~/genomics/EquSeq/mapping/wgs_mover.sh`
13. re-generate the bam list 
    - `cd $EPHEMERAL/ancestry/; ls $EPHEMERAL/wgs_data/final/*.bam > bam.list; echo $EPHEMERAL/novel_data/merged/final.bam >> bam.list`

**misc**

* Bowtie2 aligner trialled but poor results were observed


## Reference genome
1. EquCab3
2. Download from NCBI, header positions needed correcting from ncbi notation to chr:00000 format. `updateHeader.sh`
3. indexed using `bwa index`


--------
# Ancestry 
## 1. Data preparation
### 1.1 Variant calling (light)

* A light variant calling was conducted, `samtools mpileup`, done so to generate a rough list of SNPs, and marking their positions for subsequent analysis, effectively no true filtering, but sped up the angsd GL analysis
* Filtering was not harsh as this was controlled for later by `angsd`
* `samtools` and `bcftools` were used
* All sequences were used to do this 

*gatk HaplotypeCaller* - trialed but not good for multiple samples and type of calling (i.e. too high quality and too long) 

To speed up the process snps-lists were generated (rather than running on a whole chromosome)
1. get lengths of chromosomes from EquCab3 assembly report
  * local: `sh scripts/snp_caller_lists.sh`
2. generate the snp files with a set number of sites (n=10 million)
  * local: `python3 scripts/snp_caller_lists.py`
3. generate vcf files `snps.sh` 
  * `samtools mpileup -uf $REF -b $EPHEMERAL/ancestry/bam.list -l $SNPLIST | bcftools call -mv`
4. generate niave snp lists (to be run against on angsd)
  * run `ancestry/snpHandler.py -c snps` from `vcfsnps.sh` 
    - manipulates vcf files and creates a full list of niave snps
    - set the number of files to create (n = 400)
      - will correspond to the number of angsd runs (1 snp file per angsd run)

### 1.2. Genotype liklihoods (angsd) - Initial pass
1. Quality control tests conducted using `angsd` and results helped determine filtering values
2. run `ancestry.sh` 
  * `angsd` - genotype liklihoods calculated with the following:
    * `-doMaf 1`: Frequency (fixed major and minor)
    * `-GL 1` : genotype likelihood model - Samtools
    * `-doGlf 2` : beagle likelihood file dump file 
    * `-doMajorMinor 1` : assign the major and minor alleles from GL
    * `-minMaf 0.02`: Allele frequencies below 2% were removed
liklihoods of genotype for each individual site

outputs: http://www.popgen.dk/angsd/index.php/Genotype_Likelihoods (major and minor alleles, 0=A, 1=C, 2=G, 3=T)

full code:
```
$ANGSD -bam $ANC_DIR/bam.list -ref $REF -P 7 \
      -out $ANC_DIR/all_beagles/$BASE -rf $SNP \
      -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
      -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 \
      -checkBamHeaders 0 \
      -GL 1 -doGlf 2 -doMajorMinor 1  -doMaf 1 \
      -minMaf 0.02 
```


### 1.2.1. Genotype liklihoods (angsd) - bcf generation
* I also ran `-doPost 1 -doGeno 8 -doCounts 1 -doBcf 1` with the above command to generate bcf files
* Used for some of the analysis (eg. Treemix)
* using snps in final beagle output, generate a list of snps
    - `cd $EPHEMERAL/ancestry; for f in all_beagles/*.list.beagle.gz.list.out; do cat $f; echo; done > snp.all.list`
    - `scp your_path/ephemeral/ancestry/snp.all.list data/ancestry/`
    - `python3 scripts/prep_snps_vcf.py`
* rerun angsd GL but with `doBCF` to generate a vcf file
    - `ancesty.sh`
* will be used for treemix


### 1.3. Select SNPs and manipulating beagle files
* Selection of snps was conducted to minimise the influenced of linkage disequilibrium
  - i.e. snps shouldn't be too close together
1. run `window.sh` for each beagle file
  - `python $CODE_DIR/ancestry/snpHandler.py -c window -i $FILE -o $DIR/all_beagles/$BASE.list.out -d 5`
  - `-d 5` denotes snps should be 5kb apart
This will: 
i. extracted all snps from beagle file
ii. bin snps in 5kb increments
iii. 1 snp per bin is randomly sampled
iv. 5kb snp list is now used to filter beagle file 
2. merge all beagle files into 1 (to be used for agmixture and pca)
  - run `beagleMerged.sh`
3. remove przewalski horses - some of the PCA analysis
    - `beagle_no_przewalski.sh`

## 2. PCA 

1. `pcagsd` - covariance matrix based on genotype liklihoods
  - `sh bamAncestry.sh -p ALL.merged.beagle.gz`
    - tested with 0.2 and 0.5 minmaf with little difference
    - next remove prez and rerun
      - `beagle_no_przewalski.sh`
  - scp to local:`scp your_path/ephemeral/ancestry/ALL.PCA.* results/ancestry/ALL_5kb_02maf/`
  - `-minMaf 0.02`
  - `-inbreed 1` (might not be used)
2. eigen values calculated and plotting conducted in `R`
  - run `Rscript ancestry/clusters.R` to generate a clst-like file used for labelling
  - `sh ancestry/quickplot.sh -p`
    - within `plotPCA.R` switch between all individuals, and a subset excluding przewalski horses

*if issues arise when installing create use conda myenv*

## 3. Admixture
1. `ngsadmix` - estimating individual admixture proportions, based on genotype liklihoods
  - for k=2...11 (number of breeds)
  - `-minMaf 0.02`
  - `qsub -J 2-11 pca_admix.sh`
      - pops 2 to 11
2. plot `plotAdmix.R`
  - `zip -R admix 'ALL.MIX*'`
  - `scp your_path/ephemeral/ancestry/admix.zip results/ancestry/ALL_5kb_02maf/`
  - using the same cluster file as the PCA
3. `validate_K.R` used to verify which K to use
  - generates plot (included in supplementary)

## 4. TreeMix
1.  concatinate bcf files 
    - use bcf file from **2.1.**
    - `module load bcftools/1.3.1; bcftools concat -O v -o merge.vcf *.bcf`
    - snps used will be the same as those in the admixture and pca
2. convert vcf files to treemix format
    - `trees_prep_files.sh`
        - The Rscript section causes errors when submitted as a job, instead run the first part in the terminal and sub the rest as a job
    - The script:
      1. rename sample header in vcf stripping paths and extensions
      2. converts vcf to plink
      3. updates map file for correct naming of the chromosomes (required for later analysis)
      4. generate stratified freq file with plink
      5. use `plink2treemix.py` to convert plink to a treemix format
3. run treemix
  - `trees.sh` - all snps
    - `trees_benson.sh` - snps associated to benson
      - not required in the end
  - with 0 migrations
  - Przewalski as the outgroup

4. plot locally
  - `scp your_path/ephemeral/ancestry/treemix/poporder* results/ancestry/treemix/`
  - `scp your_path/ephemeral/ancestry/treemix/results/* results/ancestry/treemix/`


## 5. F3 statistic
* submit `f34_treemix.sh` as a job
* (A; B, C)
* f3 score sig. different from 0 indicates A is a result of admixture from B and C 
* I use Z <= -2 (seems to be the rule of thumb for f3)


## 6. F4 statistic
* submit `f34_treemix.sh` as a job
* test for introgression from incomplete lineage sorting
* (A,B; C,D)
* without admixture allele freq. difference in A and B should be independent from C and D (i.e. f4 = 0)
* -ve f4 indicates gene flow between C and B or A and D
* +ve f4 indicates gene flow between A and C or B and D
* Therefore, if we use A as an outgroup, knowing that there is no admixture into C or D we can test for gene flow
* +ve f4 indicates gene flow between B and D (A as outgroup)
* -ve f4 indicates gene flow between B and C (A as outgroup)
* I use a threshold as Z <= -2 and Z >= 2 (seems to be the rule of thumb for f4)
* Use przewalski as the outgroup
* Benson as target (B)
* Remove przewalski hybrid as not to appear as C or D

## Files
`ALL.MIX.pdf` : Admixture
`NO.PRZ.PCA.pdf`: PCA (0.02)
`ALL.02.PCA.pdf` : PCA (0.02 with przewalski)
`ALL.05.PCA.pdf` : PCA (0.05 with przewalski)
`tree.pdf`: maximum liklihood tree
`tree_res.pdf`: residuels from tree


--------

# Gene to trait 
## 1. Data Preparation
* Variant-trait data from:
  - Ten years of the horse reference genome: insights into equine biology, domestication and population dynamics in the post- genome era
* only used mendelian trait mutations
* only used substitutions
  - 47/56 of the total mendelian variant-traits

## Generating Data
1. generate list of sites of interest
  - `python3 phenotype_variants.py`
  - from: "Ten years of the horse reference genome: insights into equine biology, domestication and population dynamics in the post-genome era"
    - also used omia but unsure on the validity
      - SQL database setup from their sqldump 
  - sync list to hpc
2. Genotype posterior probabilities in ANGSD
  - `sh ~/genomics/EquSeq/gene_to_trait/genotypeLiklihood.sh`
    - or sub as job
  - `-doGeno 8 -doGlf 2 -dumpCounts 2  -doDepth 1 -doCounts 1`  allows for per individual allele freq counts to be conducted
    - if 0 no reads are present - assign 0
```
$ANGSD -bam $BAM_LIST -ref $REF -P 4 -out $DIR/gl.out -rf $SNP \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 \
    -C 50 -baq 1 -minMapQ 15 -minQ 15 -checkBamHeaders 0 \
    -doMajorMinor 1 -doPost 2 -doMaf 1 -minMaf 0 -gl 1 \
    -doGeno 8 -doGlf 2 -dumpCounts 2  -doDepth 1 -doCounts 1
```
3. scp GL files to local
  - `scp your_path/ephemeral/gene_to_trait/gl.out.* results/gene_to_trait/`
4. `gunzip results/gene_to_trait/gl.out.geno.gz`

## Files
`heat.pdf` : heatmap of phenotypes, and individuals 
`trait_codes.csv` : traits and markers that correspond to graph codes

--------

# Oral Microbial Diversity  
## 1. Data Preparation
### 1.1 kraken

* custom KRAKEN database 
* RefSeq for: 
  - bacteria, archaea, virus, plant, fungi, horse, human
* The header for EquCab3 was updated, allowing for use in Kraken2

1. run `refGenomes.sh` to prepare horse and human genomes
2. run `dbDownload.sh` to download the remaining data and add the above
3. run `dbbuild.sh` to build

### 1.2 taxonomic data
`taxonkit` used to assign complete lineages to kraken reports. 

1. NCBI Taxonomy database dump, containing all known phylogenies
  - 20/05/2020 downloaded


## 2. screening
  - `screen.sh`
     - `--paired` used for paired reads
     - reports generated
  - Raw fastq files (pairs) screened across db 
  - LCA mapping for each k-mer
    - taxa classified or unclassified to a k-mer
2. check confidence intervals of the clssification
  - `kraken_stat.sh` 
  - `python3 oral_diversity/kraken_explore.py -i sandbox/kraken/all.tmp.kraken -o sandbox/kraken/stats -t 0.2`
  - 20% confidence on kmer classification using C/Q scores (see kraken manual for definition)
  - pass and fail csvs will be generated
3. combine pass fail csvs updating the counts coloumn 
  - `python3 oral_diversity/taxaLineages.py -c merge -d results/oral_diversity/stats/ -o results/oral_diversity/stat.all.pass -t pass`
  - `python3 oral_diversity/taxaLineages.py -c merge -d results/oral_diversity/stats/ -o results/oral_diversity/stat.all.fail -t fail`
4. assign taxa lineages to taxa ids in stats files
  - `python3 oral_diversity/taxaLineages.py -c taxa -p results/oral_diversity/stat.all.pass.csv -f results/oral_diversity/stat.all.fail.csv -o results/oral_diversity/ -t df`
  - creates taxa id data from ncbi taxa trees 
  - impementing `taxonkit`
    - used to produce flat taxanomic lineages

## 3. Summary stats calculation
5. explore data, plot, calculate biodiversity index
  - run `microbeDiversity.R`
  - joins taxa data with counts in flat df
  - total read counts were taken from all files 
    - i.e. read1 + read2 + read3 
    - read files might be repeats of the same regions, but might also be independent
    - by averaging (and assuming only repeats) we may be biasing the data
    - eg. kmerA sampled twice, kmerB once, kmerC once  
    - most other papers are not clear what they do here, but it seems as if they take all of the reads and total them before further analysis
  - relative abundances calculated
    - as in Lassalle 2017
    - "ratio of read counts at one specific level over the total"
  - for virus and bacteria calculate:
    - shannon index
    - effective species number


**N.B.**
Kraken probably shouldn't be used to compute the abundance of species, instead Braken should be used.
Also, Lassalle 2017: "Alpha diversities were computed using phylodiversity metrics (McCoy & Matsen, 2013)"

https://www.mdpi.com/1999-4915/11/5/435/htm
*issues running? run within an active conda env*

# Files
* `kraken_db_composition.txt`: counts of species at high taxa levels (plant, fungi, etc.)
* `df_bacteria.csv` : count df for bacteria
* `df_virus.csv` : count df for virus
* `oralDiv_bacteria.pdf` : read count bar plot bacteria
* `oralDiv_virus.pdf` : read count bar plot virus
* `summary.txt` : summary stats from diversity calculations

--------


