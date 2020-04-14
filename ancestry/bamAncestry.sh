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

#cp $HOME/genomics/code/plotQC.R $TMPDIR


#FILES=$EPHEMERAL/mapping/merged/new.rg.bam
#FILES=$EPHEMERAL/mapping/merged_v1/gatk_test/file.bam
DIR=$EPHEMERAL/sra_data/
#DIR=$EPHEMERAL/mapping/merged_v1/gatk_test/

# refernece genome
REF=$EPHEMERAL/sra_data/EquCab2.fna

#cp $DIR/plotQC.R $TMPDIR


#echo '=================================='
#echo -e "\nsamtools index\n"

#samtools faidx $REF
# generates .fai

# timer
#timer



echo '=================================='
echo -e "\nANGSD\n"

# P 4 - 4 threads -maxDepth 500
# -minMapQ 20 min quality 
# -trim 0 no trimming of reads
# -maxDepth 500 : bin together depths of equal or greater

#angsd -P 4 -i $FILES -ref $REF -out $DIR/ALL.qc \
#		-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
#		-trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 13 \
#		-doQsDist 1 -doDepth 1 -doCounts 1 #\
#		#-GL 1 -doGlf 4 -doMajorMinor 1 -doMaf 1 \
		#-doGeno 32 -doPost 1


# create a list of bam files 
#ls -d $PWD/*.bam > bam.filelist

angsd -nThreads 31 -bam $DIR/bam.filelist -ref $REF -out $DIR/results/data \
		-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
		-trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 13 \
		-doQsDist 1 -doDepth 1 -doCounts 1 \
		-checkBamHeaders 0 -SNP_pval 1e-3 \
		-GL 2 -doGlf 2 -doMajorMinor 1 -doMaf 1 \
		-doGeno 32 -doPost 1



# timer
timer

echo '=================================='
echo '=================================='
echo -e "\nPCA\n"
echo '=================================='
echo '=================================='

#echo '=================================='
#echo -e "\nBuild pcangsd\n"
#git clone https://github.com/Rosemeis/pcangsd.git

#cd pcangsd/
#python setup.py build_ext --inplace

echo '=================================='
echo -e "\npcangsd - covariance matrix and admixture\n"
PCANGSD=$EPHEMERAL/sra_data/pcangsd/pcangsd.py

# Estimate covariance matrix and individual admixture proportions
python $PCANGSD -beagle $DIR/results/data.beagle.gz \
		-admix -o $DIR/results/test -threads 31




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


#echo '=================================='
#echo -e "\nplotting\n"

# counts of quality scores
#less -S $DIR/ALL.qc.qs
# counts of per-sample depth
#less -S $DIR/ALL.qc.depthSample 
#wc -l $DIR/ALL.qc.depthSample # 30 $DIR/ALL.qc.depthSample

# counts of global depth
#less -S $DIR/ALL.qc.depthGlobal

#It is convenient to compute the percentiles of these distribution
#Rscript plotQC.R $DIR/ALL.qc 

#less -S ALL.qc.mafs.gz
