#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=5:mem=7gb

# remove duplicates and poor quality reads


#----- load modules ----#
echo '=================================='
echo -e "\nLoad modules\n"
module load samtools/1.3.1
module load java/jdk-8u144
module load picard/2.6.0

#----- variables ----#
PICARD=$PICARD_HOME/picard.jar
DIR=$EPHEMERAL/test_align/




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


# Summarise alignments
echo '=================================='
echo -e "\nSummarise\n"

#java -jar $PICARD CollectAlignmentSummaryMetrics \
#          R=$EPHEMERAL/screen/EquCab2.fna \
#          I=$EPHEMERAL/merged/merged_reads.bam \
#          O=$EPHEMERAL/merged/output.txt