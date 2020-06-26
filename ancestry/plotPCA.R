# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-05-05
# Last Modified: 2020-06-26

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

MAF <- "0.05"
#covIN <- "results/ancestry/eu_more_wg_5kb_05maf/ALL.PCA.cov"
#sitesIN <- "results/ancestry/eu_more_wg_5kb_05maf/ALL.PCA.sites"
covIN <- "results/ancestry/all_test/ALL.PCA.cov"
sitesIN <- "results/ancestry/all_test/ALL.PCA.sites"
out <- "results/ancestry/ALL.PCA.pdf"

###########################################
############### Wraggling #################
###########################################

# Read sites file
nSites <- nrow(read.table(sitesIN, stringsAsFact=F))
title_append <- paste("nsites=", nSites, " ", "minMaf=", MAF)

###########################################
################## Plot ###################
###########################################

#cov <- as.matrix(read.table("results/ancestry/test.cov"))
#e <- eigen(cov)
#
#
#plot(e$vectors[,1:2], col=ID$CLUSTER)
#legend("topleft",fill=1:4,levels(ID$CLUSTER))

###### run PCA script ########
	# Modified Matteo's plotPCA.R 
		# not enough colour palette options
		# and fine tune ggplot
#com <- "Rscript scripts/plotPCA.R -i results/ancestry/test.cov -c 1-2 -a results/ancestry/test.clst -o results/ancestry/test.pca.pdf"

# exec
#system(com)
cbPalette <- c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

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
PC$Pop <- factor(annot$CLUSTER)

title <- paste("PC",comp[1]," (",signif(eig$val[comp[1]], digits=3)*100,"%)"," / PC",comp[2]," (",signif(eig$val[comp[2]], digits=3)*100,"%)",sep="",collapse="")

x_axis = paste("PC",comp[1],sep="")
y_axis = paste("PC",comp[2],sep="")

### colour palette ###
n <- length(unique(PC$Pop))
col_vector <- c(brewer.pal(8, "Dark2"), brewer.pal(8, "Set2"), brewer.pal(12, "Paired"))
col_vector <- unique(col_vector)

shapes <- c(0:n)# c(0:7, 15:18)

# sort order and levels
#poplv <- levels(PC$Pop)
#levels(PC$Pop) <- c("BENSON", sort(poplv[poplv!="BENSON"]))
#factor(PC$Pop, levels = c("BENSON", sort(poplv[poplv!="BENSON"])))
#reorder(PC$Pop, new.order=c("BENSON", sort(poplv[poplv!="BENSON"])))

#PC$PC1
#PC$PC2
g <- ggplot() + 
		geom_point(data=PC, aes_string(x=x_axis, y=y_axis, color="Pop")) + 
		ggtitle(paste(title, title_append)) + 
		theme_bw() +
		guides(col = guide_legend(ncol = 2))
		#scale_colour_manual(values = col_vector)+
		#scale_shape_manual(values =shapes) + 
		#theme(legend.position=c(1,1), 
		#	legend.justification=c(1,1))


#g <- ggplot() + 
#		geom_point(data=PC, aes_string(x=x_axis, y=y_axis, shape="Pop", color="Pop")) + 
#		ggtitle(paste(title, title_append)) + 
#		theme_bw() #+
#		scale_colour_manual(values = col_vector)+
#		scale_shape_manual(values =shapes) + 
#		#theme(legend.position=c(1,1), 
#		#	legend.justification=c(1,1))

pdf(file=out, 15, 7)
print(g)
invisible(dev.off())





