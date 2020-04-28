# specify reference and ancestral sequences
	# these can be the same for my purposes


# import unix functions
source unix_functions.sh

#echo '=================================='
#echo -e "\nLoading modules\n"
#module load angsd/915
ANGSD=$EPHEMERAL/dependencies/angsd/angsd
module load samtools/1.3.1 
module load anaconda3/personal

#echo '=================================='
#echo -e "\nDirectories\n"

DIR=$EPHEMERAL/all_data/
ANC_DIR=$EPHEMERAL/ancestry/


# refernece genome
REF=$EPHEMERAL/mapping/ref_genome/EquCab2.fna
SAMPLE=$EPHEMERAL/mapping/merged/new.rg.bam

 

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
	$ANGSD -nThreads 31 -bam $ANC_DIR/bam.list -ref $REF -out $ANC_DIR/ALL.qc\
	        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -r $CHR \
	        -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 \
	        -doQsDist 1 -doDepth 1 -doCounts 1 -checkBamHeaders 0



	echo '=================================='
	echo -e "\nplotting\n"
	cp $HOME/genomics/code/plotQC.R $TMPDIR

	#It is convenient to compute the percentiles of these distribution
	Rscript plotQC.R $ANC_DIR/ALL.qc


}

function genotypeLH(){

	CHR=$1

	# filter with the above and:
		# -doMaf 1: Frequency (fixed major and minor)
		# -doPost 1: estimate per-site allele frequency as a 
			#prior for genotype proportions
		# -GL 1 : genotype likelihood model - Samtools
		# -doGlf 2 : beagle likelihood file dump file
		# -doMajorMinor 1 : assign the major and minor alleles from GL
		# -doGeno 32 : genotype probabilities in binary format

	echo '=================================='
	echo -e "\nGenerating Genotype Liklihoods\n"

	$ANGSD -nThreads 31 -bam $ANC_DIR/bam.list -ref $REF \
			-out $ANC_DIR/ALL_$CHR -r $CHR \
			-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
			-trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 \
			-checkBamHeaders 0 \
			-GL 1 -doGlf 4


}


function genoCalling(){


	echo '=================================='
	echo -e "\nGenotype Calling\n"

	$ANGSD -glf $ANC_DIR/ALL_XXXXX.glf.gz -fai $REF.fai \
			-nInd 10 -out $ANC_DIR/ALL_CALL_$CHR \
			-doMajorMinor 1 -doGeno 3 -doPost 2 -doMaf 1


}


function pcaGL(){


	echo '=================================='
	echo -e "\npcangsd - covariance matrix and admixture\n"
	PCANGSD=$EPHEMERAL/dependencies/pcangsd/pcangsd.py

	# Estimate covariance matrix and individual admixture proportions
	python $PCANGSD -beagle $ANC_DIR/ALL.beagle.gz \
			-admix -o $ANC_DIR/ -threads 31

}


#############################
############ Main ###########
#############################

while getopts "mq:g:c:p" opt; do
  case ${opt} in
  	m) # make bam list
		makeBamList
		;;
    q ) # quality control 
       	qualityCheck $OPTARG
      	;;
    g ) # genotype liklihoods
      	genotypeLH $OPTARG
      	;;
    c) # genotype calling
		genoCalling $OPTARG
		;;
    p) # pcangsd
    	pcaGL
    	;;
    \? ) echo "Usage: cmd [-m] [-q] chrN [-g] chrN [-p]"
      ;;
  esac
done



