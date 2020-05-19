


echo '=================================='
echo -e "\nLoading modules\n"

#module load java/jdk-8u144 # gatk associated
module load anaconda3/personal
#module load gatk/4.0
module load samtools/1.3.1
module load bcftools/1.3.1
#module load picard/2.6.0 # cleaning
#PICARD=$PICARD_HOME/picard.jar

echo '=================================='
echo -e "\nDirectories\n"


DIR=$EPHEMERAL/snp_calling/
REF=$EPHEMERAL/mapping/ref_genome/EquCab2.fna
#DIC_OUT=$EPHEMERAL/mapping/ref_genome/EquCab2.dict
#SOURCE=$EPHEMERAL/mapping/merged/new.rg.bam


#function makeDic(){
#	
#	echo '=================================='
#	echo -e "\npicard reference dict\n"
#	
#	
#	#prepare the reference genome
#	java -Xmx60g -jar $PICARD CreateSequenceDictionary \
#	      R=$REF \
#	      O=$DIC_OUT\
#	      TMP_DIR=$TMP_DIR # resolves memory issues
#	
#	# timer
#	timer
#
#}

#function validate(){
#	java -jar picard.jar ValidateSamFile \
#      I=input.bam \
#      MODE=SUMMARY
#}


#function gatkCaller(){
#	# one at a time - no good 
#	### REMOVE ME ###
#
#	CHR=$1
#	FILE=SRR515214
#
#	echo '=================================='
#	echo -e "\ngatk\n"
#
#	source activate myenv # activate conda environment
#
#	# generate the VCF
#	gatk --java-options "-Xmx120g" HaplotypeCaller \
#		--reference $REF \
#		--input $EPHEMERAL/all_data/sorted/$FILE.sorted.bam \
#		--output $DIR/snps_$CHR.$FILE.vcf \
#		--intervals $CHR \
#		--tmp-dir $TMPDIR
#	timer
#
#}


function samtoolsCaller(){

	CHR=$1
	echo $CHR
	echo "out: " $DIR/snps.$CHR.raw.vcf
	echo '=================================='
	echo -e "\nsamtools\n"

	 samtools mpileup -uf $REF -b $DIR/bam.list -r $CHR \
		| bcftools call -mv > $DIR/snps.$CHR.raw.vcf


}


samtoolsCaller $1


#############################
############ Main ###########
#############################

#while getopts ":g:s" opt; do
#  case ${opt} in
#  	g) # gatk snp caller with chr
#		gatkCaller $OPTARG
#		;;
#	s) # gatk snp caller with chr
#		samtoolsCaller $OPTARG
#		;;
#    \? ) echo "Usage: cmd [-g] [-s]"
#      ;;
#  esac
#done



