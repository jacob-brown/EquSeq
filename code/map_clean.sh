#!/bin/bash
#PBS -l walltime=05:00:00
#PBS -l select=1:ncpus=32:mem=124gb


###########################################
# timer function

time_start=$SECONDS

function timer {
	duration=$(($SECONDS - $time_start))
	echo -e "\n..........................\n"
 	echo "Time elapsed: " $duration " sec"
 	echo -e "\n..........................\n"
}



###########################################
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
echo -e "\nFixmate info\n"


java -Xmx120g -jar $PICARD FixMateInformation \
		INPUT=$DIR/merged_reads.bam \
		OUTPUT=$DIR/merged_reads.fix.bam  \
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
		INPUT=$DIR/merged_reads.fix.bam \
		OUTPUT=$DIR/merged_reads.rmdup.bam  \
		M=metrics \
		REMOVE_DUPLICATES=true \
		TMP_DIR=$TMPDIR


# timer
timer

echo '-----------------------'
echo -e "\nCheck for duplicates\n"
# duplicates


samtools view -f 1024 $DIR/merged_reads.rmdup.bam | \
		wc -l > $DIR/dup_n.txt


echo "duplicates"
cat $DIR/dup_n.txt


# timer
timer

echo '-----------------------'
echo -e "\nNew file\n"

# -h to include header., -q 10 quality above 10, 
# Make a new bamfile, where you only the reads where both ends maps, 
# and filter out those with a mapping quality below 10, and removing duplicates
samtools view -h -f 2 -F 1024 $DIR/merged_reads.rmdup.bam -q 10 > $DIR/new.bam

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