# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-05-05
# Last Modified: 2020-05-12

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

covIN <- "results/ancestry/tmp_5kb/ALL.PCA.cov"
out <- "results/ancestry/ALL.PCA.pdf"

###########################################
############### Wraggling #################
###########################################



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
#levels(PC$Pop) <- c("Source", sort(poplv[poplv!="Source"]))
#factor(PC$Pop, levels = c("Source", sort(poplv[poplv!="Source"])))
#reorder(PC$Pop, new.order=c("Source", sort(poplv[poplv!="Source"])))


g <- ggplot() + 
		geom_point(data=PC, aes_string(x=x_axis, y=y_axis, shape="Pop", color="Pop")) + 
		ggtitle(title) + 
		scale_colour_manual(values = col_vector)+
		#scale_colour_discrete(breaks = leg_lvs) +
		scale_shape_manual(values =shapes) + 
		theme_bw() #+
		#theme(legend.position=c(1,1), 
		#	legend.justification=c(1,1))

pdf(file=out, 7, 5)
print(g)
invisible(dev.off())





