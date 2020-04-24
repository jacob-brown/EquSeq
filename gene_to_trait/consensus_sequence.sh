#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=32:mem=124gb

echo '=================================='
echo -e "\nLoading modules\n"
module load samtools/1.3.1 # general
module load bcftools/1.3.1


REF_GEN=$EPHEMERAL/mapping/ref_genome/EquCab2.fna
BAM=$EPHEMERAL/mapping/merged/new.rg.bam
OUT_FQ=$EPHEMERAL/gene_to_trait/cns.fq

samtools mpileup -uf $REF_GEN $BAM | \
	bcftools call -c - | \
	vcfutils.pl vcf2fq > $OUT_FQ


#samtools mpileup -su -l snp.list -f $REF_GEN new.rg.bam | bcftools view - > test.vcf

samtools mpileup -suO -l snp.list -f $REF_GEN new.rg.bam | bcftools view - > test.vcf

#samtools mpileup -s -l snp.list -f $REF_GEN new.rg.bam | \
#	bcftools view - > test1.vcf

