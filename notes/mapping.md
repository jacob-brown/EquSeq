#Mapping
 
## Tools
fastx - trim ends of reads
BWA - mapping reads onto reference genome
samtools - converting sam to bam and merging
gunzip - unzipping 

## General notes
BWA to map onto reference genome
	BWA-MEM algoritm - fast, accurate, longer reads


## Steps
1 ) get the refernce seq
	rsync from source
2) Trim ends of sequence - those with poor kmer scores
	fastx last 10 bases of all reads

---- BWA ---
3) Index  reference index 
	once, as all reads can use once indexed
4) Align to reference genome
	Remember that my Illumina data is pair ended eg. 518_1 and 518_2
	Align each pair of reads to the ref genome, generating a SAM file for each pair

	don't align them all at once, possibly fragment length penalties

---- SAMTOOLS ---
5) Convert to sam format and sort file
sorting brings the header to the top and allows for merging

6) Merge SAM files into one
`samtools merge finalBamFile.bam *.bam`


