#!/bin/bash
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=32:mem=62gb

echo '=================================='
echo -e "\nLoading modules\n"
module load fastx/0.0.14 # trimming 
module load bwa/0.7.8 # alignment
module load samtools/1.3.1 # general
module load java/jdk-8u144 # picard associated
module load picard/2.6.0 # cleaning

echo '=================================='
echo -e "\nDefine paths and files\n"


PICARD=$PICARD_HOME/picard.jar
DIR=$EPHEMERAL/mapping/
REF=$EPHEMERAL/ref_genome/EquCab2.fna

# get pairs of reads

DATA_ALL=($DIR/reads/*) # array of all data _1 and _2 extensions
FILE=${DATA[$PBS_ARRAY_INDEX]} # select the data by the job number
NOEXT=$(echo $f | cut -f 1 -d '.') # remove extension
BASE=$(basename "$NOEXT") # remove path

# index reference genome
# trim
# align
# convert to bam
# sort
# index 
# stats
# cleaning
	# fixmate
	# remove duplicates
		# check
	# mapped
	# unmapped
	# merge

echo '=================================='
echo -e "\nMove scripts to TMPDIR\n"
mv $HOME/genomics/code/map_* $TMPDIR

echo '==================================================='
echo '==================================================='
echo -e '\n---------- Begin Running Scripts ----------\n'
echo '==================================================='
echo '==================================================='

# trimming reads
	# file in and file out required
bash map_trim.sh $FILE $DIR/trimmed/$BASE.fq


echo '=================================='
echo -e "\nFixmate info\n"


java -Xmx60g -jar $PICARD FixMateInformation \
			INPUT=$DIR/read_out.sorted.bam \
			OUTPUT=$DIR/read_out.fix.bam  \
			SORT_ORDER=coordinate \
			TMP_DIR=$TMPDIR # resolves memory issues

echo '=================================='
echo -e "\nRemove Duplicates\n"

# mark duplicates
java -Xmx60g \
			-jar $PICARD MarkDuplicates \
			INPUT=$DIR/read_out.fix.bam \
			OUTPUT=$DIR/read_out.fix.md.bam  \
			M=metrics \
			REMOVE_DUPLICATES=true \
			TMP_DIR=$TMPDIR


echo '=================================='
echo -e "\nCheck for duplicates\n"
# duplicates

samtools view -f 1024 $DIR/read_out.fix.md.bam | wc -l > $DIR/duplicate_count.txt

echo '=================================='
echo -e "\nUnmapped\n"
# unmapped 
samtools view -F 4 $DIR/read_out.fix.md.bam | wc -l > $DIR/unmapped_count.txt
samtools view -F 4 $DIR/read_out.fix.md.bam > $DIR/unmapped.bam

echo '=================================='
echo -e "\nMapped\n"
# mapped 
samtools view -f 4 $DIR/read_out.fix.md.bam | wc -l > $DIR/mapped_count.txt
samtools view -f 4 $DIR/read_out.fix.md.bam > $DIR/mapped.bam







