
# vcf conversion to treemix input
# run treemix and plot results

#cd sandbox


CODEDIR=~/genomics/EquSeq/

#--------------- file prep ---------------#

echo '=================================='
echo "\ncopying vcf file to sandbox\n\n"
#cp snps/snps.chr3.raw.vcf ./snps.vcf
mv snps/snps.79.raw.vcf ./snps.vcf
#bcftools view snps/out.bcf  > snps.vcf

#cp snps/snps.vcf .
#tail snps.vcf

# KEEP FOR POSSIBLE USE LATER #
echo '=================================='
echo "\nmaking snp list\n"
#use only sites present in pca and admixture
#zcat < $CODEDIR/data/processed_sequences/beagle/ALL.merged.beagle.gz | grep chr3_ | cut -f1 | tr "_" "\t" > snp.list


echo '=================================='
echo "\ncorrecting sample names\n"
# generate a cluster file from the vcf
# get linename for and match first instance
LINE_NUM=($(grep -n -m1 CHROM snps.vcf | cut -d : -f1))

#grep -m1 CHROM snps.vcf 

# although messy, specifying linenumber speeds up the renaming process
sed $LINE_NUM's/fixflag/BENSON/g' snps.vcf | sed $LINE_NUM's!/rds/general/user/jb1919/ephemeral/wgs_data/final/!!g' | sed $LINE_NUM"s/.sorted.bam//g" | sed $LINE_NUM"s/.bam//g" > snps.rename.vcf 

# did it change the names correctly?
# grep -m1 CHROM snps.rename.vcf 

# make cluster file
	# columns: 1. Family ID; 2. ID; 3. cluster

echo '=================================='
echo "\nmaking cluster file\n"
Rscript $CODEDIR/ancestry/cluster_trees.R  snps.rename.vcf
# cat clusters.clst 
# below is modified of: https://speciationgenomics.github.io/Treemix/

echo '=================================='
echo "\nmaking treemix file\n"

# snps don't appear in map ped with plink, first use vcftools 
	# and apply some filtering
file=snps.rename

# make map and ped
	# allows for reformatting

vcftools --vcf $file".vcf" --plink --mac 2 --remove-indels --max-alleles 2 --out $file #--positions snp.list

# snps cause issues if no name is present - rename them
awk -F"\t" '{
        split($2,chr,":")
	$1="1"
	$2="1:"chr[2]
        print $0
}' ${file}.map > better.map
mv better.map ${file}.map

# generate stratified freq file 
	# and filter 
	# --mind 0.1 # important for SE but benson is absent, instead remove sites manually later
plink/plink --file snps.rename --freq --missing --snps-only --geno 0.1  --maf 0.02 --within clusters.clst --out $file #--mind 0.1 
gzip -f $file".frq.strat"
#zcat < $file".frq.strat.gz" | head


# convert using treemix python script
python2 $CODEDIR/dependancies/plink2treemix.py $file".frq.strat.gz" treemix.frq.gz
#zcat < treemix.frq.gz | head
#zcat < treemix.frq.gz | grep BENSON
#zcat < treemix.frq.gz | head -n1 | tr " " "\n" | wc -l

### remove SNPs that aren't present with our sample 
	# issues arise with treemix algorithms if lots of our sample is absent

#### KEEP FOR NOW!!!!
#python3 $CODEDIR/ancestry/bensonSNPs.py

echo "snp count:"
zcat < treemix.frq.gz | wc -l

#--------------- Analysis ---------------#
echo '=================================='
echo "\nrunning treemix\n"
treemix -i treemix.frq.gz -o tree.out.no -noss -n_warn 5 -root Przewalski > tree.out.no.log

# multiple migrations
for i in {0..5}
do 
	echo $i "/5 migrations"
 treemix -i treemix.frq.gz -m $i -o tree.out.$i -root KoninklijkWarmbloedPaardNederland -noss -n_warn 5 > tree.out.${i}.log &
done; 


# bootstrap and migration events
for i in {0..2}
do 
	echo $i "/5 migrations"
 treemix -i treemix.frq.gz -m $i -o tree.out.$i -root KoninklijkWarmbloedPaardNederland -bootstrap -k 500 -noss -n_warn 5 > tree.out.${i}_log
done; 


### repeats of the analysis ###
#for k in {0..20}
#do
#	for i in {0..5}
#	do 
#		echo $i "/5 migrations"
#	 treemix -i treemix.frq.gz -m $i -o out_treemix/tree.out.m$i$k -root SwissWarmblood -bootstrap -k 100 -noss -n_warn 5 > out_treemix/treemix_${i}_log &
#	done; wait
#done; wait


echo '=================================='
echo "\nrun threepop (treemix f3 stat)\n"
# blocks of 10 snps, as in sheep paper 
threepop -i treemix.frq.gz -k 10 > f3stat.txt
#head f3stat.txt

echo '=================================='
echo "\nplotting\n"

# create poporder file for residual plot
cat clusters.clst | cut -f3 | sort | uniq > poporder # if benson is present
cat clusters.clst | cut -f3  | sed 's/BENSON//g' | sed '/^$/d' | sort | uniq > poporder
#cat poporder

#Rscript $CODEDIR/ancestry/trees_plot.R s # single population 
#Rscript $CODEDIR/ancestry/trees_plot.R m # multiple migrations 
#Rscript $CODEDIR/ancestry/trees_plot.R l  # migration liklihood	
#Rscript $CODEDIR/ancestry/trees_plot.R f # f3 stats



echo '=================================='
echo "\ndone\n"



