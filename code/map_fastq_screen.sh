#!/bin/bash
#PBS -l walltime=00:10:00
#PBS -l select=1:ncpus=3:mem=6gb

# 0-11 done
# 12-20 doing

DATA=($EPHEMERAL/reads/*) # array of all data
FILE=${DATA[$PBS_ARRAY_INDEX]} # select the data by the job number


echo '=================================='
echo -e "\nLoad modules\n"


# do it on the trimmed data
$EPHEMERAL/screen/fastq_screen_v0.14.0/fastq_screen $FILE


echo '=================================='
echo -e "\nMove files\n"

mv *.txt $EPHEMERAL/screen/reports/
mv *.html $EPHEMERAL/screen/reports/

# gather ref genomes - remember to samtools index each
# HUNAN
# wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_cds_from_genomic.fna.gz


# E.Coli
# wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Escherichia_coli/reference/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz

# apple 
# wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/Malus_domestica/latest_assembly_versions/GCF_002114115.1_ASM211411v1/GCF_002114115.1_ASM211411v1_genomic.fna.gz

# horse2 

# wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Equus_caballus/latest_assembly_versions/GCF_002863925.1_EquCab3.0/GCF_002863925.1_EquCab3.0_genomic.fna.gz








# put config file
# scp fastq_screen.conf jb1919@login.cx1.hpc.ic.ac.uk:/rds/general/user/jb1919/ephemeral/screen/fastq_screen_v0.14.0

# grab html report
# scp jb1919@login.cx1.hpc.ic.ac.uk:/rds/general/user/jb1919/ephemeral/screen/reports/* .




