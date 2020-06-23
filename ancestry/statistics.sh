#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-10
# Last Modified: 2020-06-23
# Desc: 



cd $EPHEMERAL/sandbox
#cp ../ancestry/ALL.merged.beagle.gz .

ANGSD=$EPHEMERAL/dependencies/angsd/angsd
REALSFS=$EPHEMERAL/dependencies/angsd/misc/realSFS

REF=../ref_genome/EquCab3.fna
ANC_DIR=$EPHEMERAL/ancestry/

#### 
# let us assume there are 2 populations 

head -n23 $ANC_DIR/bam.list2 > pop1.list
tail -22 $ANC_DIR/bam.list2 > pop2.list

for i in {1..2}
do
	$ANGSD -bam pop$i.list -ref $REF -P 7 -out pop$i.out -r chr3:77000000-77005000 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -checkBamHeaders 0 -GL 1 -doGlf 1 -doMajorMinor 1 -doMaf 1 -minMaf 0.02 -doGeno 8 -doPost 1 -anc $REF -dosaf 1 
done

# estimate sfs for each pop
	# sfs, -fold 1 required if ancestor is not known
for i in {1..2}
do
	$REALSFS pop$i.out.saf.idx -fold 1 -P 4 2> /dev/null > pop$i.sfs &
done;wait

cat pop1.sfs 
cat pop2.sfs


# between populations 2dSFS
$REALSFS pop1.out.saf.idx pop2.out.saf.idx -fold 1 -P 4 2> /dev/null > pop.1.2.sfs &

cat pop.1.2.sfs

# local scp 
scp jb1919@login.cx1.hpc.ic.ac.uk:/rds/general/user/jb1919/ephemeral/sandbox/*.sfs .

# in R
R
sfs1<-scan("pop1.sfs")
sfs2<-scan("pop2.sfs")
pdf("sfs.pdf", 10, 5)
par(mfrow=c(1,2))
barplot(sfs1[-1])
barplot(sfs2[-1])
dev.off()
system("open -a Skim.app sfs.pdf")

# Fst - back on hpc
	# 2d-sfs (pop.1.2.sfs) used as a prior for the analysis
$REALSFS fst index pop1.out.saf.idx pop2.out.saf.idx -sfs pop.1.2.sfs -fstout pop.pbs -whichFST 1 &> /dev/null

# look
$REALSFS fst print pop.pbs.fst.idx

# sliding window
	# N.B. j.b. made the win and step small just for testing!!!!
$REALSFS fst stats2 pop.pbs.fst.idx -win 50 -step 10 -whichFST 1 > pop.pbs.txt

# per site fst
head pop.pbs.txt




