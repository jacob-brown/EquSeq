#! /bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=32:mem=124gb

# import unix functions

cp $HOME/genomics/code/unix_functions.sh $TMPDIR

source unix_functions.sh


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
DIR_OUT=$EPHEMERAL/gene_to_trait/
REF=$EPHEMERAL/mapping/ref_genome/EquCab2.fna
DIC_OUT=$EPHEMERAL/mapping/ref_genome/EquCab2.dict
FILES=$EPHEMERAL/mapping/merged/new.rg.bam



#echo '=================================='
#echo -e "\npicard reference dict\n"
#
#
# prepare the reference genome
#java -Xmx60g -jar $PICARD CreateSequenceDictionary \
#      R=$REF \
#      O=$DIC_OUT\
#      TMP_DIR=$TMP_DIR # resolves memory issues
#
## timer
#timer

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
gatk --java-options "-Xmx120g" HaplotypeCaller \
	--reference $REF \
	--input $FILES \
	--output $DIR_OUT/raw_variants.vcf \
	--tmp-dir $TMPDIR 	
# timer
timer