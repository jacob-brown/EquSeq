
### local ###
#cd sandbox
# generate cluster files


echo '=================================='
echo "\ncopying vcf file to sandbox\n\n"
cp snps/snps.chr3.raw.vcf ./snps.vcf


echo '=================================='
echo "\nmaking snp list\n"
# use only sites present in pca and admixture
zcat < ../data/processed_sequences/beagle/ALL.merged.beagle.gz | grep chr3_ | cut -f1 | tr "_" "\t" > snp.list


echo '=================================='
echo "\ncorrecting sample names\n"
# generate a cluster file from the vcf
# get linename for and match first instance
LINE_NUM=($(grep -n -m1 CHROM snps.vcf | cut -d : -f1))

# although messy, specifying linenumber speeds up the renaming process
sed $LINE_NUM's/fixflag/BENSON/g' snps.vcf | sed $LINE_NUM's!/rds/general/user/jb1919/ephemeral/wgs_data//sorted/!!g' | sed $LINE_NUM"s/.sorted.bam//g" > snps.rename.vcf 

# make cluster file
	# columns: 1. Family ID; 2. ID; 3. cluster

echo '=================================='
echo "\nmaking cluster file\n"
Rscript ../ancestry/cluster_trees.R  snps.rename.vcf

# below is modified from: https://speciationgenomics.github.io/Treemix/


echo '=================================='
echo "\nmaking treemix file\n"

# snps don't appear in map ped with plink, first use vcftools 
	# and apply some filtering
file=snps.rename

# make map and ped
	# allows for reformatting
#vcftools --vcf $file".vcf" --plink

vcftools --vcf $file".vcf" --plink --mac 2 --remove-indels --max-alleles 2 --out $file --positions snp.list

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
	# --mind 0.1 # important for SE but benson is absent
plink/plink --file snps.rename --freq --missing --snps-only --geno 0.1  --maf 0.02 --within clusters.clst --out $file --mind 0.1 # --mind 0.96
gzip -f $file".frq.strat"
#zcat < $file".frq.strat.gz" | head


# convert using treemix python script
python2 ../dependancies/plink2treemix.py $file".frq.strat.gz" treemix.frq.gz
#zcat < treemix.frq.gz | head
#zcat < treemix.frq.gz | grep BENSON
echo "snp count:"
zcat < treemix.frq.gz | wc -l

echo '=================================='
echo "\nrunning treemix\n"
treemix -i treemix.frq.gz -o tree.out -noss -n_warn 5 -root SwissWarmblood

# bootstrap and migration events
#for i in {0..5}
#do
# treemix -i treemix.frq.gz -m $i -o tree.out.$i -root SwissWarmblood -bootstrap -k 500 -noss -n_warn 5 > treemix_${i}_log &
#done



echo '=================================='
echo "\nplotting\n"
# create poporder file for residual plot
cat clusters.clst | cut -f3 | sort | uniq > poporder
#tail -n+2 poporder > poporder # removes benson
Rscript ../ancestry/trees_plot.R


echo '=================================='
echo "\ndone\n"



