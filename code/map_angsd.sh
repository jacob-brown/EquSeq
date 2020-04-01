#!/bin/bash
#PBS -l walltime=03:00:00
#PBS -l select=1:ncpus=5:mem=5gb


# specify reference and ancestral sequences
	# these can be the same for my purposes

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
module load angsd/915
module load samtools/1.3.1 
module load anaconda3/personal

echo '=================================='
echo -e "\nDirectories\n"

cp $HOME/genomics/code/plotQC.R $TMPDIR

#FILES=$EPHEMERAL/mapping/merged/new.rg.bam
FILES=$EPHEMERAL/sra_data/ERR868003.bam
#RES_DIR=$EPHEMERAL/mapping/merged/
RES_DIR=$EPHEMERAL/sra_data/

# refernece genome
REF=$EPHEMERAL/mapping/ref_genome/EquCab2.fna



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

angsd -i $FILES -ref $REF -out $RES_DIR/ALL.qc \
		-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
		-trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 13 \
		-doQsDist 1 -doDepth 1 -doCounts 1 #\
		#-GL 1 -doGlf 4 -doMajorMinor 1 -doMaf 1 \
		#-doGeno 32 -doPost 1

#-SNP_pval 1e-3 \
	# SNPs uncertanity,

ls $RES_DIR/ALL.*

# timer
timer

#echo '=================================='
#echo -e "\nspecific sites\n"

# add chromosome and location
#echo chr3 79504108 79618886 > $RES_DIR/snp.txt # KIT

#echo chr16 21548000 21757591 >> $RES_DIR/snp.txt # MITF

#https://www.ensembl.org/Equus_caballus/Gene/Summary?db=core;g=ENSECAG00000005360;r=16:21548000-21757591

#angsd sites index $RES_DIR/snp.txt

#angsd -i $FILES -ref $REF -out $RES_DIR/restrict.sites \
#		-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
#		-trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 13 \
#		-sites $RES_DIR/snp.txt \
#		-doQsDist 1 -doDepth 1 -doCounts 1 \
#		-GL 1 -doGlf 4 -doMajorMinor 1 -doMaf 1 \
#		-doGeno 32 -doPost 1


#echo '=================================='
#echo -e "\nSNP Calling significance\n"

# iterate over some cutoffs (you can change these)
#for PV in 0.05 1e-2 1e-4 1e-6
#do
#	if [ $PV == 0.05 ]; then echo SNP_pval NR_SNPs; fi
#	angsd -glf $RES_DIR/ALL.qc.glf.gz -nInd 1 -fai $REF.fai \
#			-out $RES_DIR/glf.$PV -doMajorMinor 1 -doMaf 1 \
#			 -skipTriallelic 1 -SNP_pval $PV &> /dev/null
#	echo $PV `zcat $RES_DIR/glf.$PV.mafs.gz | tail -n+2 | wc -l`
#done


#echo '=================================='
#echo -e "\nplotting\n"

# counts of quality scores
#less -S $RES_DIR/ALL.qc.qs
# counts of per-sample depth
#less -S $RES_DIR/ALL.qc.depthSample 
#wc -l $RES_DIR/ALL.qc.depthSample # 30 $RES_DIR/ALL.qc.depthSample

# counts of global depth
#less -S $RES_DIR/ALL.qc.depthGlobal

#It is convenient to compute the percentiles of these distribution
Rscript plotQC.R $RES_DIR/ALL.qc 

#less -S ALL.qc.mafs.gz
