#! /bin/bash
#PBS -l walltime=00:10:00
#PBS -l select=1:ncpus=5:mem=5gb
# Desc: prepare treemix files on the hpc

module load anaconda3/personal
module load vcftools
module load plink

CODEDIR=~/genomics/EquSeq/
DIR=$EPHEMERAL/ancestry/treemix/
FILE=$DIR/merge.vcf
FILE_NEW=$DIR/snps.rename

echo '=================================='
echo "\ncorrecting sample names\n"
	# generate a cluster file from the vcf
	# get linename for and match first instance
LINE_NUM=($(grep -n -m1 CHROM $FILE | cut -d : -f1))

# although messy, specifying linenumber speeds up the renaming process
sed $LINE_NUM's/fixflag/BENSON/g' $FILE | \
	sed $LINE_NUM's!/rds/general/user/jb1919/ephemeral/wgs_data/final/!!g' | \
	sed $LINE_NUM"s/.sorted.bam//g" | sed $LINE_NUM"s/.bam//g" | \
	sed $LINE_NUM's!/rds/general/user/jb1919/ephemeral/novel_data/merged/!!g' \
	> $FILE_NEW".vcf" 

# did it change the names correctly?
# grep -m1 CHROM snps.rename.vcf 

# make cluster file
	# columns: 1. Family ID; 2. ID; 3. cluster


echo '=================================='
echo "\nmaking cluster file\n"
	# issues here run manually
#Rscript $CODEDIR/ancestry/cluster_trees.R  \
#	$FILE_NEW".vcf" \
#	$CODEDIR/data/cleaned_data/info_all.csv \
#	$DIR/ALL.clst



################################################
### ALL OF THE 	ABOVE CAN BE RUN IN THE TERMINAL, SUB THE REST BELOW AS A JOB
################################################
# below is modified of: https://speciationgenomics.github.io/Treemix/


echo '=================================='
echo "\nmaking treemix file\n"

# snps don't appear in map ped with plink, first use vcftools 
	# and apply some filtering

# make map and ped
	# allows for reformatting
echo "vcftools"

vcftools --vcf $FILE_NEW".vcf" --plink --mac 2 \
	--remove-indels --max-alleles 2 --out $FILE_NEW #--positions snp.list

echo "correcting map"
# snps cause issues if no name is present - rename them
awk -F"\t" '{
        split($2,chr,":")
	$1="1"
	$2="1:"chr[2]
        print $0
}' ${FILE_NEW}.map > better.map
mv better.map ${FILE_NEW}.map

echo "generate stratified freq file"
# generate stratified freq file 
	# and filter 
	# --mind 0.1 # important for SE but benson is absent, instead remove sites manually later
plink --file $FILE_NEW --snps-only --geno 0.1 \
	--maf 0.02 --freq --missing \
	--within $DIR/ALL.clst --out $FILE_NEW #--mind 0.1 

echo "zipping"
gzip -f $FILE_NEW".frq.strat"

#zcat < $file".frq.strat.gz" | head

echo "plink2treemix.py"
# convert using treemix python script
python2 $EPHEMERAL/dependencies/plink2treemix.py $FILE_NEW".frq.strat.gz" $DIR/treemix.frq.gz
#zcat treemix.frq.gz | head
#zcat treemix.frq.gz | grep BENSON
#zcat treemix.frq.gz | head -n1 | tr " " "\n" | wc -l

### remove SNPs that aren't present with our sample 
	# issues arise with treemix algorithms if lots of our sample is absent

#### KEEP FOR NOW - remove snps not associated with BENSON
	# run as parallel analysis
python $CODEDIR/ancestry/bensonSNPs.py

echo "snp count:"
zcat $DIR/treemix.frq.gz | wc -l

echo '=================================='
echo "\nprepare cluster file for plotting\n"
# create poporder file for residual plot
cat $DIR/ALL.clst | cut -f3 | sort | uniq > poporder # if benson is present
cat $DIR/ALL.clst | cut -f3  | sed 's/BENSON//g' | sed '/^$/d' | sort | uniq > poporder.benson





