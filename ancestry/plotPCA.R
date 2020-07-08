# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-05-05
# Last Modified: 2020-07-07

# Desc: generate a table of names for running the pcangsd script


###########################################
################# Modules #################
###########################################

require(tools)
require(dplyr)
require(ggplot2)
require(RColorBrewer)

###########################################
############## Function(s) ################
###########################################



###########################################
######### Input(s) and Parameters #########
###########################################

przewalski <- T

MAF <- "0.02"
covIN <- "results/ancestry/ALL_5kb_02maf/ALL.PCA.cov"
sitesIN <- "results/ancestry/ALL_5kb_02maf/ALL.PCA.sites"
out <- "results/ancestry/ALL.PCA.pdf"

## below is for dataset with no prw horses
if(przewalski){
	covIN <- "results/ancestry/NO_PREZ_5kb_02maf/NO.PRZ.PCA.cov"
	sitesIN <- "results/ancestry/NO_PREZ_5kb_02maf/NO.PRZ.PCA.sites"
	out <- "results/ancestry/NO.PRZ.PCA.pdf"
}
###########################################
############### Wraggling #################
###########################################

# Read sites file
nSites <- nrow(read.table(sitesIN, stringsAsFact=F))
title_append <- paste("nsites=", nSites, " ", "minMaf=", MAF)

###########################################
################## Plot ###################
###########################################


###### run PCA script ########
	# Modified Matteo's plotPCA.R 
		# not enough colour palette options
		# and fine tune ggplot
#com <- "Rscript scripts/plotPCA.R -i results/ancestry/test.cov -c 1-2 -a results/ancestry/test.clst -o results/ancestry/test.pca.pdf"

# exec
#system(com)
cbPalette <- c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbPalettelrg <- c("#004949","#009292","#ff6db6","#ffb6db",
 					"#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
 					"#920000","#924900","#db6d00","#24ff24","#ffff6d", "#000000")
# Read input file
covar <- read.table(covIN, stringsAsFact=F);

# Read annot file
annot <- read.table("results/ancestry/clusters", sep="\t", header=T) # note that plink cluster files are usually tab-separated

# Parse components to analyze
comp <- as.numeric(strsplit("1-2", "-", fixed=TRUE)[[1]])

# Eigenvalues
eig <- eigen(covar, symm=TRUE);
eig$val <- eig$val/sum(eig$val);
cat(signif(eig$val, digits=3)*100,"\n");

# Plot
PC <- as.data.frame(eig$vectors)
colnames(PC) <- gsub("V", "PC", colnames(PC))

#####
# remove przewalski? - testing a different PCA 
if(covIN == "results/ancestry/NO_PREZ_5kb_02maf/NO.PRZ.PCA.cov"){
	annot <- annot[(annot$CLUSTER != "Przewalski" & annot$CLUSTER != "Przewalski-hybrid"),]
}

# update cluster values
PC$Pop <- factor(annot$CLUSTER)

title <- paste("PC",comp[1]," (",signif(eig$val[comp[1]], digits=3)*100,"%)"," / PC",comp[2]," (",signif(eig$val[comp[2]], digits=3)*100,"%)",sep="",collapse="")

x_axis = paste("PC",comp[1],sep="")
y_axis = paste("PC",comp[2],sep="")

popCount <- length(unique(PC$Pop))

# ggrepel::geom_text_repel() + 
# theme(legend.position = "none") 
#geom_text(position=position_jitter(width=1,height=1)) + 
g <- ggplot(data=PC, aes_string(x=x_axis, y=y_axis, color="Pop", label = "Pop")) + 
		geom_text() + 
		ggtitle(paste(title, title_append)) + 
		theme_bw() +
		scale_colour_manual(values = rep(cbPalette, len = popCount))+
		theme(legend.position = "none") 
		#guides(col = guide_legend(ncol = 2))

pdf(file=paste0(out,".txt.pdf"), 15, 15)
print(g)
invisible(dev.off())

p <- ggplot(data=PC, aes_string(x=x_axis, y=y_axis, color="Pop", shape = "Pop")) + 
		geom_point(size = 3) + 
		ggtitle(paste(title, title_append)) + 
		theme_bw() +
		scale_shape_manual(values = rep(0:18, len = popCount)) + 
		scale_colour_manual(values = rep(cbPalette, len = popCount))+
		guides(col = guide_legend(ncol = 2))+
		theme(legend.position="bottom")

pdf(file=out, 15, 15)
print(p)
invisible(dev.off())


system(paste0("open -a Skim.app ", out,".txt.pdf"))
system(paste0("open -a Skim.app ", out))



