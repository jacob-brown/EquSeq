#!/bin/bash
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=1:mem=10gb

# path to files
PATH_IN='/rds/general/user/jb1919/home/genomics/sequences/cdts-hk.genomics.cn/Clean/F19FTSEUHT1854-swab-horse-1A/'

# reference genome
#REF_GEN='/rds/general/user/jb1919/home/genomics/sequences/refgenome/EquCab2.fna'

#F1='/rds/general/user/jb1919/home/genomics/sequences/cdts-hk.genomics.cn/Clean/F19FTSEUHT1854-swab-horse-1A/V300044309_L2_B5GHORlfyRAAAAAAA-517_1.fq'

#F2='/rds/general/user/jb1919/home/genomics/sequences/cdts-hk.genomics.cn/Clean/F19FTSEUHT1854-swab-horse-1A/V300044309_L2_B5GHORlfyRAAAAAAA-517_2.fq'

# path to output
#FILE_OUT='/rds/general/user/jb1919/home/genomics/sequences/mapped/ali.sam'

# list of files
DATA=$(ls $PATH_IN/*.fq)

#----- copy to temp ----#
# reads - 85GB of data
#echo '=================================='
#echo -e "\nmoving reads to temp\n"
#rsync -rltv $DATADIR $EPHEMERAL
#rsync -rltv $F1 $TMPDIR/read_1.fq
#rsync -rltv $F2 $TMPDIR/read_2.fq

#echo '=================================='
#echo -e "\nmoving ref genome to temp\n"

# ref genome
#rsync -rltv $HOME/genomics/sequences/refgenome/EquCab2.fna $TMPDIR/EquCab2.fna 

echo '=================================='
echo -e "\nmoving indexed ref genome to temp\n"
rsync -rltv /rds/general/user/jb1919/home/genomics/sequences/refgenome/indexed/* $TMPDIR/ 

#----- Load Modules ----#
# load BWA 
echo '=================================='
echo -e "\nLoading BWA\n"
module load bwa/0.7.8

#----- index ----#
# index reference genome
#echo '=================================='
#echo -e "\nIndex ref genome\n"
#bwa index EquCab2.fna
# roughly 30 mins

# move indexed ref genome and rename
#echo '=================================='
#echo -e "\nmoving indexed ref gen\n"
#mv EquCab2.* $HOME/genomics/sequences/refgenome/indexed/

#----- align ----#

# align sequences
echo '=================================='
echo -e "\nAlign sequences\n"
#bwa mem ref.fa read1.fq read2.fq > aln-pe.sam
bwa mem EquCab2.fna $DATA > align.sam
#bwa mem EquCab2.fna read_1.fq.gz read_2.fq.gz > align.sam


echo '=================================='
echo -e "\nMove sam sequences\n"
mv align.sam $HOME/genomics/sequences/mapped/align.sam






