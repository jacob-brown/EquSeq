


echo '=================================='
echo -e "\nLoading modules\n"

module load anaconda3/personal
module load samtools/1.3.1
module load bcftools/1.3.1


echo '=================================='
echo -e "\nDirectories\n"


DIR=$EPHEMERAL/snp_calling/
#REF=$EPHEMERAL/ref_genome/EquCab3.fna
REF=~/genomics/old_wd/ref_genome/EquCab3.fna



function samtoolsCaller(){

	SNPLIST=$1
	JOB_N=$2
	OUT=$DIR/snps.$JOB_N.raw.vcf
	echo $SNPLIST
	echo "out: " $OUT
	echo '=================================='
	echo -e "\nsamtools\n"


	samtools mpileup -uf $REF -b $EPHEMERAL/ancestry/bam.list -l $SNPLIST \
		| bcftools call -mv > $OUT

}


samtoolsCaller $1 $2


# old methods - 1 chrom per run
#	 samtools mpileup -uf $REF -b $EPHEMERAL/ancestry/bam.list -r $CHR \
#		| bcftools call -mv > $DIR/snps.$CHR.raw.vcf

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



