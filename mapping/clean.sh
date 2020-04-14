#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=32:mem=124gb


# import unix functions
source ../scripts/unix_functions.sh


echo '=================================='
echo -e "\nLoading modules\n"

module load samtools/1.3.1 # general
module load java/jdk-8u144 # picard associated
module load picard/2.6.0 # cleaning

# remove duplicates and poor quality reads

DIR=$EPHEMERAL/mapping/merged/
PICARD=$PICARD_HOME/picard.jar

# timer
timer

echo '-----------------------'
echo -e "\nCheck for duplicates\n"
# duplicates

echo "duplicates: "
samtools view -f 1024 --threads 31 $DIR/merged.bam | wc -l 


echo '-----------------------'
echo -e "\nFixmate info\n"
# check mate-pair information is in sync between each read and its mate pair

java -Xmx120g -jar $PICARD FixMateInformation \
		INPUT=$DIR/merged.bam \
		OUTPUT=$DIR/merged.fix.bam  \
		SORT_ORDER=coordinate \
		TMP_DIR=$TMPDIR # resolves memory issues


# timer
timer
#
echo '-----------------------'
echo -e "\nRemove Duplicates\n"

# mark duplicates
java -Xmx120g \
		-jar $PICARD MarkDuplicates \
		INPUT=$DIR/merged.fix.bam \
		OUTPUT=$DIR/merged.rmdup.bam  \
		M=metrics \
		REMOVE_DUPLICATES=true \
		TMP_DIR=$TMPDIR


# timer
timer

echo '-----------------------'
echo -e "\nCheck for duplicates\n"
# duplicates

echo "duplicates: "
samtools view -f 1024 --threads 31 $DIR/merged.rmdup.bam | wc -l 


# timer
timer

echo '-----------------------'
echo -e "\nNew file\n"

# -h to include header., -q 20 quality above 20, 
# Make a new bamfile, where you only the reads where both ends maps, 
# and filter out those with a mapping quality below 20, and removing duplicates
samtools view -h -f 2 -F 1024 --threads 31 \
		$DIR/merged.rmdup.bam -q 20 > \
		$DIR/new.bam

# timer
timer


echo '-----------------------'
echo -e "\nAdding RG flag\n"

# fix the read groups
java -Xmx120g -jar $PICARD AddOrReplaceReadGroups \
	I=$DIR/new.bam \
	O=$DIR/new.rg.bam \
	RGLB=lib1 \
	RGPL=illumina \
	RGPU=unit1 \
	RGSM=fixflag \
	SORT_ORDER=coordinate \
	CREATE_INDEX=True \
	TMP_DIR=$TMPDIR # resolves memory issues

echo '-----------------------'
echo -e "\nIndex with samtools\n"

samtools index $DIR/new.rg.bam

# timer
timer