# EquSeq

`git clone ...`


Always run the scripts from the parent directory `EquSeq`, relative paths have been used.

## Languages
* python
* R
* bash

## Software and Downloads
* fastqc - quality scores
* FASTX - trimming the last 10 bases
* BWA - aligning and mapping
* picard - marking reads for duplicates and updating info
* samtools - general read cleaning and processing
* kraken2 - metagenomic analysis of fastq reads
* taxonkit
* angsd
* pcagsd

### Conda environment

## Files
* `info_all.csv` - Unique sequence runs/file codes with additional information. 
* `info_individual_grouped.csv` - As above, but grouped at the individual sample level. Noting the number of files on NCBI associated with the BioSample ID.
* `info_pop_grouped.csv` - As above, but grouped at the sub-group level. Noting how many individual samples we have for each species and sub_group.

Collection - runinfo tables manually from run selector pasting in all bioproject codes
* `expand_projects.py` - expand the project codes and prepare list for a deeper search.


### Data

* `SraRunTable.txt` lists all the metadata from the following NCBI query: `(((horse[Organism]) OR Equus[Organism]) AND genomic[Source]) AND WGS[Strategy]`
* `SRR_Acc_List.txt` is a list of the `Accession Codes`
* `prj_modern.csv` - bioproject codes for modern sequences
* `prj_ancient.csv` - bioproject codes for ancient sequences


### File structure

run: `tree -I 'sandbox' -d > tree.md`

```
.
├── ancestry
├── data
│   ├── ancestry
│   │   ├── snp.chr
│   │   └── test_success_beagle
│   ├── cleaned_data
│   │   └── infotables_update
│   ├── gene_variants
│   ├── metadata
│   ├── omia_sql
│   ├── oral_diversity
│   │   └── taxonkit_misc
│   ├── processed_sequences
│   │   ├── sandbox_data
│   │   └── snps
│   └── raw_data
│       ├── infotables_original
│       └── supplementary_data_from_studies
├── gene_to_trait
├── job_submissions
├── mapping
├── oral_diversity
├── results
│   ├── ancestry
│   │   ├── wg_5kb_02maf
│   │   └── wg_5kb_05maf
│   └── oral_diversity
│       └── kraken_reports
├── scripts
└── setup

```