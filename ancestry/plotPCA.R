# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-05-05
# Last Modified: 2020-07-24

# Desc: generate a table of names for running the pcangsd script


###########################################
################# Modules #################
###########################################

require(tools)
require(dplyr)
require(ggplot2)
require(RColorBrewer)
require(ggforce)
require(gghighlight)
require(ggpubr)
require(scales)

###########################################
############## Function(s) ################
###########################################

readPC <- function(covIN, sitesIN, MAF, clustIn, breedIN){


	# Read sites file
	nSites <- nrow(read.table(sitesIN, stringsAsFact=F))
	title_append <- paste("nsites=", nSites, " ", "minMaf=", MAF)
	title_append <- paste("minMaf=", MAF)

	# Read input file
	covar <- read.table(covIN, stringsAsFact=F)

	# Read annot file
	annotAll <- read.table(clustIn, sep="\t", header=T) 

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
		annot <- annotAll[(annotAll$CLUSTER != "Przewalski" & annotAll$CLUSTER != "Przewalski-hybrid"),]
	}else{
		annot <- annotAll
	}

	# levels
	clst <- sort(as.character(unique(annotAll$CLUSTER)))
	prw <- c("Przewalski-hybrid", "Przewalski")
	clust_levels <- c(clst[!(clst %in% prw)], prw)

	# update cluster values
	PC$Pop <- factor(annot$CLUSTER, levels=clust_levels)
	PC$IID <- as.character(annot$IID)
	title <- paste("PC1"," (",signif(eig$val[1], digits=3)*100,"%)"," / PC",
				"PC2"," (",signif(eig$val[2], digits=3)*100,"%)",sep="",collapse="")
	title <- paste0(title, " ", title_append)
	# abbreviate
	PC$Abv <- as.vector(sapply(as.character(PC$Pop), function(x) substr(x, 1, 3)))

	# add major categories
	breeds <- read.csv(breedIN, header=F)
	breeds$V2[breeds$V2 == ""] <- NA
	colnames(breeds) <- c("Pop", "maj_grp")

	PC <- left_join(PC, breeds, by="Pop")

	# assign shape (for consistant legend when prw are removed)
	#popCount <- length(unique(annot$CLUSTER))
	#shapeVal <- c(15:18, 3, 4, 8 ,10)
	#shapes <- rep(shapeVal, len = popCount)
	#shapedf <- data.frame(shape=shapes, Pop=unique(annot$CLUSTER))
	#PC <- left_join(PC, shapedf, by="Pop")

	#issueIDs <- c("32", "90", "120")
	issueIDs <- c("90")
	PC <- PC[!PC$IID %in% issueIDs,] # issue ID

	return(list(PC, title))
}

plotPCA <- function(PC_ls, zoom=T, highlight=NA, arrow=F, 
	leg_pos="right", n_cols=1, nbreak=3, title=NA){

	PC <- PC_ls[[1]]
	if(is.na(title)){
		title <- PC_ls[[2]]
	}

	# counts
	popCount <- length(unique(PC$Pop))
	
	# colours
	cbPalette <- c("#000000", "#999999", "#E69F00", "#56B4E9", 
				"#009E73", "#F0E442", 
				"#0072B2", "#D55E00", "#CC79A7")
	# shapes
	shapes <- c(15:18, 3, 4, 8 ,10)

	p <- ggplot(data=PC, aes(x=PC1, y=PC2, 
				color=Pop, shape = Pop, fill=Pop)) + 
			geom_point(size = 3,  alpha = 0.85, stroke = 0.7) + 
			ggtitle(title) + 
			scale_shape_manual(values = rep(shapes, len = popCount)) + 
			scale_colour_manual(values = rep(cbPalette, len = popCount))+
			scale_x_continuous(breaks=trans_breaks(identity, identity, n = nbreak))+
			theme_pubr()+
			labs(shape = "Breeds", color = "Breeds", fill = "Breeds") +
			theme(legend.position=leg_pos,
					legend.key.size = unit(3, "mm"),
					panel.grid.major = element_blank(),
					panel.grid.minor = element_blank(),
					panel.background = element_blank()
					)+
			guides(shape=guide_legend(ncol=n_cols,bycol=TRUE,
							title.position="top")) 							

	if(zoom){
		p <- p + 
				facet_zoom(ylim=c(-0.04,0.12), zoom.size=2.5, horizontal = T)
	}

	if(arrow){
		p <- p +
			annotate(
    			geom = "curve", 
    			#x = xs, y = ys, 
    			#xend = xe, yend = ye, 
    			xend = -0.0039, yend = 0.025, 
    			x = -0.015, y = 0.015, 
    			curvature = 0, arrow = arrow(length = unit(2, "mm"))
  				) 
	}

	if(!is.na(highlight)){
		p <- p + gghighlight(Pop %in% highlight)
	}

	return(p)

}

###########################################
######### Input(s) and Parameters #########
###########################################

clustPath <- "results/ancestry/clusters"
breedPath <- "results/ancestry/breeds.csv"

###########################################
############### Wraggling #################
###########################################

# 0.02 minmaf - No przewalski
nop05 <- readPC(covIN = "results/ancestry/NO_PREZ_5kb_02maf/NO.PRZ.PCA.cov", 
				sitesIN = "results/ancestry/NO_PREZ_5kb_02maf/NO.PRZ.PCA.sites", 
				MAF = "0.02", 
				clustIn = clustPath,
				breedIN = breedPath)

### misc plots ###
# All individuals 0.02 minmaf
all02 <- readPC(covIN = "results/ancestry/ALL_5kb_02maf/ALL.02.PCA.cov", 
				sitesIN = "results/ancestry/ALL_5kb_02maf/ALL.02.PCA.sites", 
				MAF = "0.02", 
				clustIn = clustPath,
				breedIN = breedPath)

# All individuals 0.05 minmaf
all05 <- readPC(covIN = "results/ancestry/ALL_5kb_05maf/ALL.05.PCA.cov", 
				sitesIN = "results/ancestry/ALL_5kb_05maf/ALL.05.PCA.sites", 
				MAF = "0.05", 
				clustIn = clustPath,
				breedIN = breedPath)

###########################################
################## Plot ###################
###########################################

### plot
outpca <- "results/ancestry/pca.pdf"
outleg <- "results/ancestry/pca_leg.pdf"

# 1
p1 <- plotPCA(PC_ls=all02, zoom=F, 
				leg_pos="bottom", 
				n_cols=3, nbreak=5, arrow=F, title="A")
# 2
p2 <- plotPCA(PC_ls=nop05, zoom=T, leg_pos="none", 
				n_cols=3, arrow=F, title="B")
# arrange plots
arrangep <- ggarrange(p1, p2, ncol=1, nrow=2, legend="none")

# get legend and plot as seperate plot
leg <- get_legend(p1)

# save 
ggsave(filename=outpca, 
	plot=arrangep, device ="pdf", width=160, height=200,  units="mm", dpi ="screen")
ggsave(filename=outleg, 
	plot=as_ggplot(leg), device ="pdf", 
	width=160, height=60,  units="mm", dpi ="screen")
system(paste0("open -a Skim.app ", outpca))
system(paste0("open -a Skim.app ", outleg))



### minmaf 05 (supplementary)
p3 <- plotPCA(PC_ls=all05,  zoom=F, leg_pos="bottom", n_cols=3)
out3 <- "results/ancestry/ALL_05_PCA.pdf"
ggsave(filename=out3, 
	plot=p3, device ="pdf",width=180, height=180, units="mm", dpi ="screen")
system(paste0("open -a Skim.app ", out3))









