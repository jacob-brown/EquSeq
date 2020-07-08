#! /bin/bash
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=1:mem=1gb


#cp -r $EPHEMERAL/oral_diversity/db_kraken_horse/ $EPHEMERAL/oral_diversity/bkup_db_kraken_horse



#mkdir novel_merged/ snp_calling dependencies


#cp -r $EPHEMERAL/ref_genome/* ~/genomics/old_wd/ref_genome/

#mv -r $EPHEMERAL/novel_data/* # 649G
#cp $EPHEMERAL/novel_data/merged/final.bam* ~/genomics/old_wd/novel_merged/
#mv -r $EPHEMERAL/snp_calling/* #18G
#cp $EPHEMERAL/snp_calling/{snps.chr26.raw.vcf,snps.chr3.raw.vcf,snps.chr4.raw.vcf} ~/genomics/old_wd/snp_calling/
#mv -r $EPHEMERAL/eu_wgs_data/sorted/ #1.5T

# to sort
#cp -r $EPHEMERAL/dependencies/ dependencies/


# allow to go
#/rds/general/user/jb1919/ephemeral/eu_wgs_data/
#/rds/general/user/jb1919/ephemeral/oral_diversity/



# qsub -J 0-171 copy.sh
	# qsub -J 0-3 copy.sh
	# qsub -J 4-171 copy.sh

DIR=/rds/general/user/jb1919/ephemeral/wgs_data/final/
FILES=($(ls -v $DIR/*.bam))
FILE="${FILES[$PBS_ARRAY_INDEX]}"

echo "moving: " $FILE*

cp $FILE* /rds/general/user/jb1919/projects/savolainen-archive-2018/live/rawdata/Horses/bam_files/



#24GB in 10mins
#1 file ~15 mins