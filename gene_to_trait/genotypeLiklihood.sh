#! /bin/bash
#PBS -l walltime=00:30:00
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
BAM_LIST=$EPHEMERAL/old_ancestry/bam.list2

# refernece genome
#REF=$EPHEMERAL/ref_genome/EquCab3.fna
REF=~/genomics/old_wd/ref_genome/EquCab3.fna

echo '=================================='
echo -e "\nGLs\n"

#$ANGSD -bam $BAM_LIST -ref $REF -P 4 \
#		-out $DIR/chr3.gl.out -r chr3:79504300-79593715 \
#		-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
#		-trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 \
#		-checkBamHeaders 0 \
#		-GL 1 -doGlf 2 -doMajorMinor 1  -doMaf 1 \
#		-minMaf 0


SNP=$CODE_DIR/data/gene_variants/trait.snps/trait.snp.mend.list
### with a bcf output ###

# beagle
$ANGSD -bam $BAM_LIST -ref $REF -P 4 -out $DIR/gl.out -rf $SNP -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -checkBamHeaders 0 -doMajorMinor 1 -doPost 1 -doMaf 1 -minMaf 0 -gl 1 -doGeno 8 -doGlf 2


# bcf
#$ANGSD -bam $BAM_LIST -ref $REF -P 4 -out $DIR/gl.out -rf $SNP -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 15 -minQ 15 -checkBamHeaders 0 -doMajorMinor 1 -doPost 1 -doMaf 1 -minMaf 0 -gl  1 -dobcf 1 --ignore-RG 0 -docounts 1 -doGeno 8 

#$ANGSD -beagle gl.out.beagle.gz -fai $REF.fai -out gl.new -doGeno 4


# chr3:79593650-79593670
# chr16:20103081-20105348


#$ANGSD

# local
# scp 


##### local testing #####
#dependancies/angsd/angsd -bam data/processed_sequences/benson/final.bam -ref data/processed_sequences/ref_genome/EquCab3.fna -P 4 -out sandbox/gl.out -r data/gene_variants/trait.snps/trait.snp.9.list -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 15 -minQ 15 -checkBamHeaders 0 -doMajorMinor 1 -doPost 1 -doMaf 1 -minMaf 0 -gl 1 -dobcf 1 --ignore-RG 0 -doGeno 1 -docounts 1






