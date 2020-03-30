#!/bin/bash
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=5:mem=15gb

module load tabix
module load igvtools/2.3.57

IGV=$IGVTOOLS_HOME/bin/igvtools.jar

java -Xmx12g $IGVTOOLS_HOME/bin/igvtools.jar index $EPHEMERAL/mapping/old_merged/raw_variants_RG.vcf

java -Xmx12g $IGV sort -t tmp $EPHEMERAL/mapping/old_merged/raw_variants_RG.vcf \
				$EPHEMERAL/mapping/old_merged/var.sort.vcf

gzip $EPHEMERAL/mapping/old_merged/var.sort.vcf

#tabix -p $EPHEMERAL/mapping/old_merged/raw_variants_RG.vcf.gz

