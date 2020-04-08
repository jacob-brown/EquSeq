#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-04-07
# Last Modified: 2020-04-07
# Desc: 

# ERR2179543.bam - Welsh Pony 
# ERR868003.bam - Shetland
# file.bam - our sample

# write a table for IDing the sequences
	# ensure order is the same as the bamlist?
	
# Rscript -e should work...
#'write.table(cbind(seq(1,3),rep(1,3),c(rep("WEL",1),rep("SHT",1),rep("VIN",1))), row.names=F, sep="\t", col.names=c("FID","IID","CLUSTER"), file="../results/ancestry/test.clst", quote=F)'

Rscript plotPCA.R -i ../results/ancestry/test.cov -c 1-2 \
		-a ../results/ancestry/test.clst \
		-o ../results/ancestry/test.pca.pdf
