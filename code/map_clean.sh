#!/bin/bash
# remove duplicates and poor quality reads

DIR=$EPHEMERAL/mapping/
PICARD=$PICARD_HOME/picard.jar

# catch input files 
BASE_NAME=$1 

echo '-----------------------'
echo -e "\nFixmate info\n"


java -Xmx60g -jar $PICARD FixMateInformation \
		INPUT=$DIR/sorted/$BASE_NAME'.sorted.bam' \
		OUTPUT=$DIR/tmp_cleaning/$BASE_NAME'.fix.bam'  \
		SORT_ORDER=coordinate \
		TMP_DIR=$TMPDIR # resolves memory issues

echo '-----------------------'
echo -e "\nRemove Duplicates\n"

# mark duplicates
java -Xmx60g \
		-jar $PICARD MarkDuplicates \
		INPUT=$DIR/tmp_cleaning/$BASE_NAME'.fix.bam' \
		OUTPUT=$DIR/tmp_cleaning/$BASE_NAME'.fix.md.bam'  \
		M=metrics \
		REMOVE_DUPLICATES=true \
		TMP_DIR=$TMPDIR


echo '-----------------------'
echo -e "\nCheck for duplicates\n"
# duplicates

samtools view -f 1024 $DIR/tmp_cleaning/$BASE_NAME'.fix.md.bam' | \
		wc -l > $DIR/stats/$BASE_NAME'.dup_n.txt'

echo "duplicates"
cat $DIR/stats/$BASE_NAME'.dup_n.txt'

echo '-----------------------'
echo -e "\nUnmapped\n"

# unmapped 
samtools view -F 4 $DIR/tmp_cleaning/$BASE_NAME'.fix.md.bam' > \
		$DIR/tmp_cleaning/$BASE_NAME'.unmapped.bam'

echo '-----------------------'
echo -e "\nMapped\n"

# mapped 
samtools view -f 4 $DIR/tmp_cleaning/$BASE_NAME'.fix.md.bam' > \
		$DIR/cleaned/$BASE_NAME'.mapped.bam'





