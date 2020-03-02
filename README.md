# Horse Genomics

## Software and Downloads
Used
* fastqc - quality scores
* FASTX - trimming the last 10 bases
* BWA - aligning and mapping

To Look into
* Guy: blobtools,Busco software, spades

possibly used
* SRA Toolkit
* Entrez Direct


## Running To Do 
* List all datasets (bioproject codes)
* script to find data associated with bioproject codes
* script to filter any unusuable data, summary table of data types and sample size
* Can other data be gathered from the list of all bioproject codes?



## Workflow
### Overview
If BAM files are used the pipline and methods should be the same, as otherwise we won't be able to compare the sequences. 

FASTQ would allow us to use our own mapping and pipeline, however will be alot more work. 

For the writeup we should aim for a hypothesis driven approach. We use the sample and want to find out the ancestry, mixture, etc. Rather than a methods write-up.

### Steps
**1.Available Data**
* Look at Orlando
* Collect all available data (everything!)
* What sort of format is it in? FASTQ, BAM?
* What pipeline/protocols where used? are they comparable?
* Pipeline and protocol should be comparable. 
* Where is it from?

**1.1 Ancient DNA (extra)**
* needs to be anlysied seperatly from the modern, and should be done after the gathering of the modern samples. 
* Common issues: Less coverage, potential contamination, DNA damage, transition/transversion tends to be an issue, centre of fragments is usually good but edges are not. 
* When and where is it from?

**2.Pipeline/Protocol**
* from Orlando - how was the data processed?
Mapping softwares - samtools etc. 

**3. Analysis**
* AMSD, NGSTOOLS, etc. 
* blobtools, Busco software, spades, BWA
* Admixture plots, PCAs, trees, etc. 
* lots can be done. 


## Files
Collection - runinfo tables manually from run selector pasting in all bioproject codes


### Code
* `expand_projects.py` - expand the project codes and prepare list for a deeper search.

### Data

* `SraRunTable.txt` lists all the metadata from the following NCBI query: `(((horse[Organism]) OR Equus[Organism]) AND genomic[Source]) AND WGS[Strategy]`
* `SRR_Acc_List.txt` is a list of the `Accession Codes`
* `prj_modern.csv` - bioproject codes for modern sequences
* `prj_ancient.csv` - bioproject codes for ancient sequences