#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-10
# Last Modified: 2020-06-26
# Desc: 



cd $EPHEMERAL/sandbox

module load anaconda3/personal

#cp ../ancestry/ALL.merged.beagle.gz .
# cp angsd_play/out.mafs.gz .
ANGSD=$EPHEMERAL/dependencies/angsd/angsd
REALSFS=$EPHEMERAL/dependencies/angsd/misc/realSFS
#REALSFS=../dependancies/angsd/misc/realSFS

#REF=../ref_genome/EquCab3.fna
REF=~/genomics/old_wd/ref_genome/EquCab3.fna

###### REMOVE head - on full run #######
# all bam lists
ALL_LISTS=($(ls ~/genomics/EquSeq/data/ancestry/bam_list_grps/ | head))

# get the snps used in the previous angsd analysis
zcat out.mafs.gz | cut -f1-2 | tail -n+2 > snp.stat.list
head snp.stat.list > tmp.list ### REMOVE LATER
wc -l snp.stat.list
wc -l tmp.list
cat tmp.list
# determine which file paths correspond to pops
# locally run: python3 ../ancestry/which_pop.py 
for list in "${ALL_LISTS[@]}"
do
	full_path=~/genomics/EquSeq/data/ancestry/bam_list_grps/$list
	pop=($(echo $list | cut -f1 -d'.' ))
	out=angsd_stats/$pop.out
	
	echo -e "\nin list: " $full_path
	echo -e "out file: "  $out "\n"
	
	# run angsd - with doSaf
	$ANGSD -bam $full_path -ref $REF -P 7 -out $out -r chr3:77000000-77005000 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -checkBamHeaders 0 -GL 1 -doGlf 1 -doMajorMinor 1 -doMaf 1 -minMaf 0.02 -doGeno 8 -doPost 1 -anc $REF -dosaf 1 > angsd_stats/log.$pop;

	# run realSFS 
	$REALSFS $out.saf.idx -fold 1 -P 4 2> /dev/null > angsd_stats/$pop.sfs
done


# local scp 
cd sandbox/stats
scp jb1919@login.cx1.hpc.ic.ac.uk:/rds/general/user/jb1919/ephemeral/sandbox/angsd_stats/*.sfs .


# plot SFSs
Rscript ../../ancestry/statistics_plot.R 

# between populations 2dSFS
# j.b. - write a script to do this for all combinations and save the results


python ~/genomics/EquSeq/ancestry/fst.py 






######### ALL re-written in python #########
pop1=akhalteke
pop2=arabian 
pop3=connemara
cd angsd_stats/
# between populations 2dSFS
# pop1 and pop2 
$REALSFS $pop1.out.saf.idx $pop2.out.saf.idx -fold 1 -P 4 2> /dev/null > $pop1.$pop2.sfs 
# pop1 and pop3 
$REALSFS $pop1.out.saf.idx $pop3.out.saf.idx -fold 1 -P 4 2> /dev/null > $pop1.$pop3.sfs 
# pop2 and pop3
#$REALSFS $pop3.out.saf.idx $pop2.out.saf.idx -fold 1 -P 4 2> /dev/null > $pop3.$pop2.sfs 

cat $pop1.$pop2.sfs 
cat $pop1.$pop3.sfs
#cat $pop3.$pop2.sfs 

# calculate Fst
	# this can be improved with multiple pops in one  !!!!!!!
$REALSFS fst index $pop1.out.saf.idx $pop2.out.saf.idx -sfs $pop1.$pop2.sfs -fstout $pop2.pbs -whichFST 1 &> /dev/null

# look at results
$REALSFS fst print $pop2.pbs.fst.idx 

# sliding window
		# N.B. j.b. made the win and step small just for testing!!!!
#$REALSFS fst stats2 $pop2.pbs.fst.idx  -win 50 -step 10 -whichFST 1 > $pop2.pbs.txt

# look at per site
cat $pop2.pbs.txt

# summary stats for the pair
$REALSFS fst stats $pop2.pbs.fst.idx > $pop2.fst.res
cat $pop2.fst.res # FST.Unweight and Fst.Weight



######### ALL re-written in python ######### To here....






















################## TESTING ################ 
##### 
## let us assume there are 2 populations 
#
#head -n23 $ANC_DIR/bam.list2 > pop1.list
#tail -22 $ANC_DIR/bam.list2 > pop2.list
#
##$ANGSD -bam $ANC_DIR/bam.list2 -ref $REF -P 7 -out pop$i.out -r chr3:77000000-77005000 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -checkBamHeaders 0 -GL 1 -doGlf 1 -doMajorMinor 1 -doMaf 1 -minMaf 0.02 -doGeno 8 -doPost 1 -anc $REF
#
#
#for i in {1..2}
#do
#	$ANGSD -bam pop$i.list -ref $REF -P 7 -out pop$i.out -r chr3:77000000-77005000 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -checkBamHeaders 0 -GL 1 -doGlf 1 -doMajorMinor 1 -doMaf 1 -minMaf 0.02 -doGeno 8 -doPost 1 -anc $REF -dosaf 1 
#done
#
#
#
## estimate sfs for each pop
#	# sfs, -fold 1 required if ancestor is not known
#for i in {1..2}
#do
#	$REALSFS pop$i.out.saf.idx -fold 1 -P 4 2> /dev/null > pop$i.sfs &
#done;wait
#
#cat pop1.sfs 
#cat pop2.sfs
#
## local scp 
#scp jb1919@login.cx1.hpc.ic.ac.uk:/rds/general/user/jb1919/ephemeral/sandbox/*.sfs .
#
## in R
#R
#sfs1<-scan("pop1.sfs")
#sfs2<-scan("pop2.sfs")
#pdf("sfs.pdf", 10, 5)
#par(mfrow=c(1,2))
#barplot(sfs1[-1])
#barplot(sfs2[-1])
#dev.off()
#system("open -a Skim.app sfs.pdf")
#
## between populations 2dSFS
#$REALSFS pop1.out.saf.idx pop2.out.saf.idx -fold 1 -P 4 2> /dev/null > pop.1.2.sfs 
#
#cat pop.1.2.sfs
#
## plot
#R
#sfs<-scan("pop.1.2.sfs")
#
#
## Fst - back on hpc
#	# 2d-sfs (pop.1.2.sfs) used as a prior for the analysis
#$REALSFS fst index pop1.out.saf.idx pop2.out.saf.idx -sfs pop.1.2.sfs -fstout pop.pbs -whichFST 1 &> /dev/null
#
## look
#$REALSFS fst print pop.pbs.fst.idx
#
## sliding window
#	# N.B. j.b. made the win and step small just for testing!!!!
#$REALSFS fst stats2 pop.pbs.fst.idx -win 50 -step 10 -whichFST 1 > pop.pbs.txt
#
## per site fst
#head pop.pbs.txt
#$REALSFS fst stats pop.pbs.fst.idx # get stats for the pair! 
#	# now use this to make a plot
#
#
# merge
#realSFS cat
