


echo '=================================='
echo -e "\nLoading modules\n"

module load anaconda3/personal
module load samtools/1.3.1
module load bcftools/1.3.1


echo '=================================='
echo -e "\nDirectories\n"


DIR=$EPHEMERAL/snp_calling/
REF=$EPHEMERAL/mapping/ref_genome/EquCab2.fna



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



