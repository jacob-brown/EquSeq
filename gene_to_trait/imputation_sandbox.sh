#! /bin/bash
#PBS -l walltime=04:00:00
#PBS -lselect=1:ncpus=32:mem=62gb

ANGSD=$EPHEMERAL/dependencies/angsd/angsd
module load anaconda3/personal
module load samtools/1.3.1

DIR=$EPHEMERAL/sandbox/


# refernece genome
REF=$DIR/EquCab3.fna

#samtools faidx $DIR/EquCab3.fna

# make bam list
#ls -d $DIR/*.bam > $DIR/bam.list

# make snp list
#snp.list.tmp

# gl

$ANGSD -bam $DIR/bam.list -ref $REF -P 7 \
			-out $DIR/out -rf $DIR/snp.list.tmp \
			-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
			-trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 \
			-checkBamHeaders 0 -doPost 1\
			-GL 1 -doMajorMinor 1  -doMaf 1 \
			-doPlink 2 -doGeno -4 \
			-minMaf 0.02

#wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20200428.zip
#unzip plink_linux_x86_64_20200428.zip

#./plink --tfile out --maf 0.02 --make-bed --out binary_fileset 