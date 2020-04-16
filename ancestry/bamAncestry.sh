#!/bin/bash
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=32:mem=62gb


# specify reference and ancestral sequences
	# these can be the same for my purposes


# import unix functions
source ../scripts/unix_functions.sh



echo '=================================='
echo -e "\nLoading modules\n"
module load angsd/915
module load samtools/1.3.1 
module load anaconda3/personal

echo '=================================='
echo -e "\nDirectories\n"

DIR=$EPHEMERAL/all_data/
ANC_DIR=$EPHEMERAL/ancestry/

# refernece genome
REF=$EPHEMERAL/mapping/ref_genome/EquCab2.fna


echo '=================================='
echo '=================================='
echo -e "\nANGSD\n"
echo '=================================='
echo '=================================='

# create a list of bam files 
#echo "Create bam list from: " $DIR/sorted/ "to: " $ANC_DIR
#ls -d $DIR/sorted/*.bam > $ANC_DIR/bam.list





function qualityCheck(){
	# check quality of bam files
		# P 4 - 4 threads -maxDepth 500
		# -minMapQ 20 min quality 
		# -trim 0 no trimming of reads
		# -maxDepth 500 : bin together depths of equal or greater
	
	echo '=================================='
	echo -e "\nChecking Quality \n"

	angsd -P 4 -b $ANC_DIR/bam.list -ref $REF -out $ANC_DIR/ALL.qc -r 11\
	        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
	        -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 \
	        -doQsDist 1 -doDepth 1 -doCounts 1 -checkBamHeaders 0



	echo '=================================='
	echo -e "\nplotting\n"
	cp $HOME/genomics/code/plotQC.R $TMPDIR

	#It is convenient to compute the percentiles of these distribution
	Rscript plotQC.R $ANC_DIR/ALL.qc


}

function genotypeLH(){

	# filter with the above and:
		# -doMaf 1: Frequency (fixed major and minor)
		# -doPost 1: estimate per-site allele frequency as a 
			#prior for genotype proportions
		# -GL 2 : genotype likelihood model - GATK
		# -doGlf 2 : beagle likelihood file dump file
		# -doMajorMinor 1 : assign the major and minor alleles from GL
		# -doGeno 32 : genotype probabilities in binary format

	echo '=================================='
	echo -e "\nGenerating Genotype Liklihoods\n"

	angsd -nThreads 31 -bam $DIR/bam.filelist -ref $REF \
			-out $ANC_DIR/ALL \
			-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
			-trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 \
			-doQsDist 1 -doDepth 1 -doCounts 1 \
			-checkBamHeaders 0 -SNP_pval 1e-3 \
			-GL 2 -doGlf 2 -doMajorMinor 1 -doMaf 1 \
			-doGeno 32 -doPost 1

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

while getopts ":qgp" opt; do
  case ${opt} in
    q ) # quality control
       	qualityCheck
      ;;
    g ) # genotype liklihoods
      genotypeLH
      ;;
    p) # pcangsd
    	pcaGL
    ;;
    \? ) echo "Usage: cmd [-q] [-g] [-p]"
      ;;
  esac
done









#echo '=================================='
#echo -e "\nspecific sites\n"

# add chromosome and location
#echo chr3 79504108 79618886 > $DIR/snp.txt # KIT

#echo chr16 21548000 21757591 >> $DIR/snp.txt # MITF

#https://www.ensembl.org/Equus_caballus/Gene/Summary?db=core;g=ENSECAG00000005360;r=16:21548000-21757591

#angsd sites index $DIR/snp.txt

#angsd -i $FILES -ref $REF -out $DIR/restrict.sites \
#		-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
#		-trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 13 \
#		-sites $DIR/snp.txt \
#		-doQsDist 1 -doDepth 1 -doCounts 1 \
#		-GL 1 -doGlf 4 -doMajorMinor 1 -doMaf 1 \
#		-doGeno 32 -doPost 1


#echo '=================================='
#echo -e "\nSNP Calling significance\n"

# iterate over some cutoffs (you can change these)
#for PV in 0.05 1e-2 1e-4 1e-6
#do
#	if [ $PV == 0.05 ]; then echo SNP_pval NR_SNPs; fi
#	angsd -glf $DIR/ALL.qc.glf.gz -nInd 1 -fai $REF.fai \
#			-out $DIR/glf.$PV -doMajorMinor 1 -doMaf 1 \
#			 -skipTriallelic 1 -SNP_pval $PV &> /dev/null
#	echo $PV `zcat $DIR/glf.$PV.mafs.gz | tail -n+2 | wc -l`
#done



