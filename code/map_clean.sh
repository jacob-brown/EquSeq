#!/bin/bash
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=5:mem=15gb

# remove duplicates and poor quality reads


#----- load modules ----#
echo '=================================='
echo -e "\nLoad modules\n"
module load anaconda3/personal
module load java
module load picard


PICARD=$PICARD_HOME/picard.jar

#samtools mpileup $EPHEMERAL/merged/merged_reads.bam | cut -f4 | sort -n | uniq -c > $EPHEMERAL/merged/dep1

#samtools rmdup -s $EPHEMERAL/merged/merged_reads.bam $EPHEMERAL/merged/merged_reads.md.bam

#samtools mpileup $EPHEMERAL/merged/merged_reads.bam | cut -f4 | sort -n | uniq -c > $EPHEMERAL/merged/dep2

echo '=================================='
echo -e "\nFix mate info\n"

# reorder and fill in the mate information
java -jar $PICARD FixMateInformation \
			INPUT=$EPHEMERAL/merged/merged_reads.bam \
			OUTPUT=$EPHEMERAL/merged/id.fixmate.srt.bam \
			SORT_ORDER=coordinate

echo '=================================='
echo -e "\nMark Duplicates\n"

# mark duplicates
java -jar $PICARD MarkDuplicates \
			INPUT=$EPHEMERAL/merged/id.fixmate.srt.bam \
			OUTPUT=$EPHEMERAL/merged/id.fixmate.srt.md.bam  \
			M=metrics


# Summarise alignments
echo '=================================='
echo -e "\nSummarise\n"

#java -jar $PICARD CollectAlignmentSummaryMetrics \
#          R=$EPHEMERAL/screen/EquCab2.fna \
#          I=$EPHEMERAL/merged/merged_reads.bam \
#          O=$EPHEMERAL/merged/output.txt