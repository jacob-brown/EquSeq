#! /bin/bash
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=5:mem=5gb


# qsub -J 0-9 genotypeLiklihood.sh


echo '=================================='
echo -e "\nLoading modules\n"

ANGSD=$EPHEMERAL/dependencies/angsd/angsd
module load anaconda3/personal

CODE_DIR=$HOME/genomics/EquSeq/
source $CODE_DIR/scripts/unix_functions.s

DIR=$EPHEMERAL/gene_to_trait/

# use same bam list as ancestry
BAM_LIST=$EPHEMERAL/ancestry/bam.list

# refernece genome
REF=$EPHEMERAL/ref_genome/EquCab3.fna

echo '=================================='
echo -e "\nGLs\n"

$ANGSD -bam $BAM_LIST -ref $REF -P 4 \
		-out $DIR/chr3.gl.out -r chr3:79504300-79593715 \
		-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
		-trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 \
		-checkBamHeaders 0 \
		-GL 1 -doGlf 2 -doMajorMinor 1  -doMaf 1 \
		-minMaf 0