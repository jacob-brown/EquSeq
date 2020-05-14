# EquSeq

`git clone `


Always run the scripts from the parent directory, i.e. from within EquSeq.

## Software and Downloads
Used
* fastqc - quality scores
* FASTX - trimming the last 10 bases
* BWA - aligning and mapping
* picard - marking reads for duplicates and updating info
* samtools - general read cleaning and processing
* kraken2 - metagenomic analysis of fastq reads

To Look into
* Guy: blobtools,Busco software, spades

possibly use
* SRA Toolkit

Conda environment
* conda install -c hcc aspera-cli



## Running To Do 
* List all datasets (bioproject codes)
* script to find data associated with bioproject codes
* script to filter any unusuable data, summary table of data types and sample size
* Can other data be gathered from the list of all bioproject codes?

## Files
* `info_all.csv` - Unique sequence runs/file codes with additional information. 
* `info_individual_grouped.csv` - As above, but grouped at the individual sample level. Noting the number of files on NCBI associated with the BioSample ID.
* `info_pop_grouped.csv` - As above, but grouped at the sub-group level. Noting how many individual samples we have for each species and sub_group.

Collection - runinfo tables manually from run selector pasting in all bioproject codes


### Code Files
* `getall_` - get all the available sequences phase of work
* `map_` - mapping phase of work



* `expand_projects.py` - expand the project codes and prepare list for a deeper search.

### Mapping 
Script run order: 
1) map_environment
2) check the reference sequence, if the IDs are not chrX then run map_correct_header.sh 
3) map_index
4) map_master
5) map_merge
6) map_clean


### Data

* `SraRunTable.txt` lists all the metadata from the following NCBI query: `(((horse[Organism]) OR Equus[Organism]) AND genomic[Source]) AND WGS[Strategy]`
* `SRR_Acc_List.txt` is a list of the `Accession Codes`
* `prj_modern.csv` - bioproject codes for modern sequences
* `prj_ancient.csv` - bioproject codes for ancient sequences


### File structure

[]: # (tree -I 'sandbox|written_project' -d -o tree.md)

```

.
├── ancestry
├── data
│   ├── cleaned_data
│   │   └── infotables_update
│   ├── gene_variants
│   ├── metadata
│   └── raw_data
│       ├── infotables_original
│       └── supplementary_data_from_studies
├── gene_to_trait
├── job_submissions
├── mapping
├── oral_diversity
├── scripts
└── setup



```