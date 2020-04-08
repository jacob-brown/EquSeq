#! /bin/bash
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=32:mem=62gb

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

#module load java/jdk-8u144 # gatk associated
#module load gatk/4.0

# remove duplicates and poor quality reads
#GATK=$GATK_HOME/bin/gatk

module load samtools/1.3.1
module load anaconda3/personal


module load java/jdk-8u144 # picard associated
module load picard/2.6.0 # cleaning
PICARD=$PICARD_HOME/picard.jar

echo '=================================='
echo -e "\nDirectories\n"

#FILES=$EPHEMERAL/mapping/merged/merged_reads.bam #merged_reads.rmdup.mapped.bam
#DIR=$EPHEMERAL/mapping/merged/
#REF=$EPHEMERAL/mapping/ref_genome/EquCab2.fna
#DIC_OUT=$EPHEMERAL/mapping/ref_genome/EquCab2.dict

DIR=$EPHEMERAL/mapping/merged/
REF=$EPHEMERAL/mapping/ref_genome/EquCab2.fna
DIC_OUT=$EPHEMERAL/mapping/ref_genome/EquCab2.dict
FILES=$EPHEMERAL/mapping/merged/new.rg.bam


#echo '=================================='
#echo -e "\nsamtools\n"

#samtools faidx $REF
#samtools sort -m 40GiB --threads 32 $FILES -o  \
		#$FILES'.sorted.bam'

#timer

echo '=================================='
echo -e "\npicard reference dict\n"


# prepare the reference genome
java -Xmx60g -jar $PICARD CreateSequenceDictionary \
      R=$REF \
      O=$DIC_OUT\
      TMP_DIR=$TMP_DIR # resolves memory issues

# timer
timer

#echo '=================================='
#echo -e "\nreordering\n"


#java -Xmx60g -jar $PICARD ReorderSam \
#	  I=$DIR/new.rg.bam \
#	  O=$DIR/reordered.bam \
#	  R=$REF 

# timer
#timer

#java -Xmx48g -jar $PICARD CleanSam \
#		I=$FILES \
#		O=$FILES'.cleaned.bam'

#echo '=================================='
#echo -e "\nIndex\n"

#samtools index  $DIR/reordered.bam

#samtools index  $FILES

#timer

echo '=================================='
echo -e "\ngatk\n"

source activate myenv # activate conda environment

# generate the VCF
#gatk HaplotypeCaller \
#	--reference $REF \
#	--input $DIR/reordered.bam \
#	--output $DIR/raw_variants.vcf \
#	--tmp-dir $TMPDIR \
#	--intervals chr3 
gatk HaplotypeCaller \
	--reference $REF \
	--input $FILES \
	--output $DIR/raw_variants.vcf \
	--tmp-dir $TMPDIR \
	--intervals chr3 	
# timer
timer