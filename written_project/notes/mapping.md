#Mapping

fastqc report suggests good base quality
BGI notes that it is phred33 

## Tools
fastx - trim ends of reads
BWA - mapping reads onto reference genome
samtools - converting sam to bam and merging
gunzip - unzipping 

## General notes
BWA to map onto reference genome
	BWA-MEM algoritm - fast, accurate, longer reads


## Steps
1 ) get the reference seq
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

7) picard: update mate data and mark duplicates

8) samtools - remove dupliactes 

9) remove unmapped reads



# Misc
Bowtie2 trialled but poor results were observed



# dump
next) remove duplicates
$SAMTOOLS rmdup -s NA19238.bam NA19238.md.bam


java -jar $PICARD MarkDuplicates I=id.fixmate.srt2.bam \
O=id.fixmate.srt.md.bam  M=metrics;

next) remove unmapped?
samtools view -h -F 4 -b blah.bam > blah_only_mapped.bam


next) remove poorly quality mapped reads

maybe with
samtools view -f 2 -F 1024 id.fixmate.srt.md.bam -q 10 >new.bam



