# specify reference and ancestral sequences
	# these can be the same for my purposes


echo '=================================='
echo -e "\nLoading modules\n"

ANGSD=$EPHEMERAL/dependencies/angsd/angsd
NGSADMIX=$EPHEMERAL/dependencies/angsd/misc/NGSadmix
PCANGSD=$EPHEMERAL/dependencies/pcangsd/pcangsd.py
module load anaconda3/personal


DIR=$EPHEMERAL/wgs_data/
ANC_DIR=$EPHEMERAL/ancestry/


# refernece genome
#REF=$EPHEMERAL/ref_genome/EquCab3.fna
REF=~/genomics/old_wd/ref_genome/EquCab3.fna
SAMPLE=$EPHEMERAL/novel_data/merged/final.bam
 

function makeBamList(){

	echo '=================================='
	echo -e "\nMaking BAM list\n"
	echo "Create bam list from: " $DIR/sorted/ "to: " $ANC_DIR
	ls -d $DIR/sorted/*.bam > $ANC_DIR/bam.list

	echo "Adding novel samples to list."
	echo $SAMPLE >> $ANC_DIR/bam.list

}	


function qualityCheck(){

	CHR=$1
	# check quality of bam files
		# P 4 - 4 threads -maxDepth 500
		# -minMapQ 20 min quality 
		# -trim 0 no trimming of reads
		# -maxDepth 500 : bin together depths of equal or greater
	
	echo '=================================='
	echo -e "\nChecking Quality \n"
	# rmove -r chr1 after
	$ANGSD -P 31 -bam $ANC_DIR/bam.list -ref $REF -out $ANC_DIR/ALL.qc\
	        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -r $CHR \
	        -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 \
	        -doQsDist 1 -doDepth 1 -doCounts 1 -checkBamHeaders 0



	echo '=================================='
	echo -e "\nplotting\n"
	cp $HOME/genomics/code/plotQC.R $TMPDIR

	#It is convenient to compute the percentiles of these distribution
	Rscript plotQC.R $ANC_DIR/ALL.qc


}


function GLBeagle(){
	# create beagle file

	SNP=$1
	BASE=($( basename $SNP))
	echo $BASE

	# filter with the above and:
		# -doMaf 1: Frequency (fixed major and minor)
		# -GL 1 : genotype likelihood model - Samtools
		# -doGlf 2 : beagle likelihood file dump file
		# -doMajorMinor 1 : assign the major and minor alleles from GL

		# -doPost 1: estimate per-site allele frequency as a 
			#prior for genotype proportions
		# -doGeno 32 : genotype probabilities in binary format

	echo '=================================='
	echo -e "\nGenerating Genotype Liklihoods\n"

	$ANGSD -bam $ANC_DIR/bam.list -ref $REF -P 7 \
			-out $ANC_DIR/all_beagles/$BASE -rf $SNP \
			-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
			-trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 \
			-checkBamHeaders 0 \
			-GL 1 -doGlf 2 -doMajorMinor 1  -doMaf 1 \
			-minMaf 0.02

}

function GLVCF(){
	# create BCF file

	SNP=$1
	BASE=($( basename $SNP))
	echo $BASE

	# doBcf 

	echo '=================================='
	echo -e "\nGenerating Genotype Liklihoods\n"

	$ANGSD -bam $ANC_DIR/bam.list -ref $REF -P 7 \
			-out $ANC_DIR/all_bcf/$BASE -rf $SNP \
			-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
			-trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 \
			-checkBamHeaders 0 \
			-GL 1 -doGlf 2 -doMajorMinor 1  -doMaf 1 \
			-minMaf 0.02 -doPost 1 -doGeno 8 -doCounts 1 -doBcf 1 
}

function pca(){


	echo '=================================='
	echo -e "\npcangsd - covariance matrix\n"
	

	BEAGLE_FILE=$1
	source activate myenv # issues with install

	# Estimate covariance matrix and individual admixture proportions
	# default -minMaf 0.05
 	python $PCANGSD -beagle $BEAGLE_FILE  \
					-o $ANC_DIR/ALL.PCA -threads 4 \
					-minMaf 0.02 -inbreed 1 -sites_save

	conda deactivate # issues with install
}


function admix(){

	echo '=================================='
	echo -e "\nNGSadmix\n"

	# default -minMaf 0.05
	BEAGLE_FILE=$1
	K_ARRAY=(2 3 4 5 6) 

	for K in "${K_ARRAY[@]}"; do

		$NGSADMIX -likes $BEAGLE_FILE -K $K\
			-outfiles $ANC_DIR/ALL.MIX.K$K -P 4 -minMaf 0.02 -seed 212
	done
	

}


#############################
############ Main ###########
#############################

while getopts "mq:g:c:p:a:b:" opt; do
  case ${opt} in
  	m) # make bam list
		makeBamList
		;;
    q ) # quality control 
       	qualityCheck $OPTARG
      	;;
    g ) # genotype liklihoods - beagle
      	GLBeagle $OPTARG
      	;;
    b) # genotype liklihoods - bcf
      	GLVCF $OPTARG
      	;;
    c) # genotype calling
		genoCalling $OPTARG
		;;
    p) # pcangsd
    	pca $OPTARG
    	;;
    a) #admixture
		admix $OPTARG
    	;;
    \? ) echo "Usage: cmd [-m] [-q] chrN [-g] chrN [-p] [-b]"
      ;;
  esac
done



