#! /bin/bash
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=5:mem=5gb


# qsub genotypeLiklihood.sh


echo '=================================='
echo -e "\nLoading modules\n"

ANGSD=$EPHEMERAL/dependencies/angsd/angsd
#module load anaconda3/personal

CODE_DIR=$HOME/genomics/EquSeq/
source $CODE_DIR/scripts/unix_functions.sh

DIR=$EPHEMERAL/gene_to_trait/

# use same bam list as ancestry
BAM_LIST=$EPHEMERAL/ancestry/bam.list

# refernece genome
#REF=$EPHEMERAL/ref_genome/EquCab3.fna
REF=~/genomics/old_wd/ref_genome/EquCab3.fna

echo '=================================='
echo -e "\nGLs\n"

SNP=$CODE_DIR/data/gene_variants/trait.snps/trait.snp.mend.all.list # combined list

# beagle, posterior probabilities from uniform dist

$ANGSD -bam $BAM_LIST -ref $REF -P 4 -out $DIR/gl.out -rf $SNP \
		-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 \
		-C 50 -baq 1 -minMapQ 15 -minQ 15 -checkBamHeaders 0 \
		-doMajorMinor 1 -doPost 2 -doMaf 1 -minMaf 0 -gl 1 \
		-doGeno 8 -doGlf 2 -dumpCounts 2  -doDepth 1 -doCounts 1

timer
#### testing 
#head $SNP > tmp.snp.list
#cat tmp.snp.list

#$ANGSD -bam $BAM_LIST -ref $REF -P 4 -out $DIR/test -rf tmp.snp.list -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 15 -minQ 15 -checkBamHeaders 0 -doMajorMinor 1 -doPost 2 -doMaf 1 -minMaf 0 -gl 1 -doGeno 8 -doGlf 2 -dumpCounts 2  -doDepth 1 -doCounts 1

#zcat test.counts.gz | cut -f172 | head 

# perhaps don't filter at angsd, rather keep the summary file and filter during the analysis stages







