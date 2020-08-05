# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-07-13
# Last Modified: 2020-08-05

# Desc: 


###########################################
################# Modules #################
###########################################

require(pophelper)
require(stringr)
require(tools)
library(gridExtra)
library(grid)
require(tidyverse)
library(ggpubr)


###########################################
############## Function(s) ################
###########################################

correctNames <- function(popls){
	
	popstore <- as.character(popls)

	correct_df <- read.csv("data/ancestry/clusters_alt.csv", header=F, stringsAsFactors=F)
	colnames(correct_df) <- c("Pop", "new")

	to_change <- correct_df[correct_df$Pop != correct_df$new,]

	runs <- length(popstore)
	for(i in 1:runs){

		if(popstore[i] %in% to_change$Pop){
			popstore[i] <- to_change[to_change$Pop == popstore[i],2] # change the name
		}
	}

	return(as.factor(popstore))
}

###########################################
######### Input(s) and Parameters #########
###########################################

maxAnc <- 9

DIR <- "results/ancestry/ALL_5kb_02maf/"
out <- "results/ancestry/admixture.pdf"
#out <- "results/ancestry/admixture_all.pdf"
clusters <- "results/ancestry/clusters"

# select the files
files <- list.files(DIR, full.names=T)
files_to_use <- subset(files, file_ext(files) == "qopt")
files_to_use <- str_sort(files_to_use, numeric = TRUE) # natural sort
files_to_use <- files_to_use[2:maxAnc]

# clusters
pop <- read.delim(clusters, header = TRUE, sep="\t", colClasses = "character")
#sapply(pop, is.character)

cbPalette <- c( "#999999", "#E69F00", "#56B4E9", 
				"#009E73", "#F0E442", "#0072B2", "#D55E00", 
				"#CC79A7", "#000000")

cbPalettelrg <- c("#004949","#009292","#ff6db6","#ffb6db",
 					"#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
 					"#920000","#924900","#db6d00","#24ff24","#ffff6d", "#000000")

#wes_pal <- wesanderson::wes_palette("Zissou1", maxAnc, type = "continuous")

###########################################
############### Wraggling #################
###########################################
# read files to qlist 
qlist_all <- readQ(files_to_use)
slist <- alignK(qlist_all, type="auto") # sort "within", "auto", "across"
k_labs <- c(paste0("K=", (1:length(slist))+1))
names(slist) <- c(paste0("K", (1:length(slist))+1))

###
# manually correct the pops that weren't autocorrected
# K3 
	# c3 > c2
names(slist$K3) <- c("Cluster1", "Cluster3", "Cluster2")
# K4 
	# c1 > c2
	# c4 > c1
names(slist$K4) <- c("Cluster2", "Cluster4", "Cluster3", "Cluster1")
# K5
	# c2 > c1
	# c1 > c5
	# c3 > c2
names(slist$K5) <- c("Cluster2", "Cluster3", "Cluster5", "Cluster4", "Cluster1")
# K6
	# c6 > c2 
	# c2 > c3
	# c5 > c1
	# c4 > c6
names(slist$K6) <- c("Cluster4", "Cluster3", "Cluster5", "Cluster6", "Cluster1", "Cluster2")
# K7
	# c3 ><c7
names(slist$K7) <- c("Cluster1", "Cluster2", "Cluster7", 
							"Cluster4", "Cluster5", "Cluster6", "Cluster3")
# K8
	# c3 ><c7
names(slist$K8) <- c("Cluster1", "Cluster2", "Cluster7", 
							"Cluster4", "Cluster5", "Cluster6", "Cluster3", "Cluster8")
# K9
	# c3 ><c7
names(slist$K9) <- c("Cluster1", "Cluster2", "Cluster7", 
							"Cluster4", "Cluster5", "Cluster6", 
							"Cluster3", "Cluster8", "Cluster9")

# clusters
onelabset1 <- pop[,3,drop=FALSE]
#onelabset1$CLUSTER <- as.vector(sapply(onelabset1$CLUSTER, function(x) substr(x, 1, 3)))
colnames(onelabset1) <- "grp"
onelabset1[onelabset1 == "BEN",] <- " BEN "
onelabset1[onelabset1 == "BENSON",] <- " BENSON "

# plot
#plotAdmix <- function(){
##showindlab=T,
#p1 <- plotQ(slist,
#			imgoutput="join",returnplot=T,exportplot=F,quiet=T,basesize=11,
#			showlegend=T, legendkeysize=10,legendtextsize=10,
#			splab=k_labs,
#			clustercol=cbPalette,
#			legendmargin=c(2,2,2,0),legendrow=1, 
#			grplab=onelabset1,#grplabsize=4,linesize=0.8,pointsize=4,
#			grplabangle=45,
#			#grplabangle=45, grplabjust=0.5, 
#			#panelratio=c(3,5), grplabheight = 10,
#			grplabsize=3.5,pointsize=6,linesize=7,linealpha=0.2,
#            pointcol="white",grplabpos=0.5,linepos=0.5,grplabheight=0.75,
#			ordergrp=T)
#
#pdf(file="sandbox/mix.pdf", 20, 10)
#grid.arrange(p1$plot[[1]],ncol=1,nrow=1)
#invisible(dev.off())
#system("open -a Skim.app sandbox/mix.pdf")
#}
#plotAdmix()

###########################################
############### Analysis ##################
###########################################

# transform to long format
store <- data.frame()
for(i in 1:length(slist)){
	tmp <- slist[[i]]
	kval <- ncol(tmp)
	tmp$X <- as.numeric(rownames(tmp))
	tmp_lng <- pivot_longer(tmp, cols= starts_with('Cluster'),
					names_to = "clst", values_to = "val") %>%
				arrange(clst, X) %>%
				mutate(K=as.factor(kval))
	store <- rbind(store, tmp_lng)
}
#head(store)
# add clusters
onelabset1$X <- as.numeric(rownames(onelabset1))
store_join <- left_join(store, onelabset1, by ="X")
store_join$X <- as.factor(store_join$X)
store_join$grp <- as.factor(store_join$grp)


# K order and update name
Klvls <- paste0("K=", as.character(unique(store_join$K)))
store_join$K <- factor(paste0("K=",as.character(store_join$K)), levels=Klvls)

#store_join$K <- as.factor(paste0("K=",as.character(store_join$K)))

# number grps 
store_join$grpN <- as.numeric(store_join$grp)

# remove cluster names
#store_join$clst
store_join$clst <- as.vector(sapply(store_join$clst, function(x) str_remove(x, "Cluster")))

# create point df to plot
pointdf <- store_join %>%
			filter(K == "K=2") %>%
			select(grp, val, X, clst) 

pointdf$val <- 0

# rename
store_join$grp <- correctNames(store_join$grp)


plotAdmix <- function(){
g <- ggplot(store_join, aes(fill=clst, colour=clst, y=val, x=X, label=clst))+
		#facet_grid(rows = vars(K))+
		facet_grid(vars(K), vars(grp), 
				scales = "free", space = "free", 
				switch="x") + 
		geom_bar(position=position_stack(reverse = TRUE), 
				stat="identity", width=1,
				size=0.01
				)+
		scale_y_continuous(expand=c(0,0))+
		#scale_y_continuous(limits = c(0,1.01), expand = expansion(mult = c(.1, .1))) +
		#theme_pubr()+
		xlab("")+
		ylab("")+
		scale_fill_manual(values = cbPalette) +
		scale_colour_manual(values = cbPalette) +
		#scale_fill_manual(values = wes_pal)+
		theme(axis.text.x = element_blank(),
				axis.line.y = element_blank(),
				axis.line.x = element_blank(),
				axis.title.y = element_blank(),
        		axis.text.y = element_blank(),
        		axis.ticks.y = element_blank(),
        		axis.ticks.x = element_blank(),
        		legend.position = "top",
        		legend.justification = "right",
				legend.key.size = unit(3, "mm"),
				strip.text.y = element_text(angle=0),
				strip.text.x = element_text(angle = 90, size = 4, hjust=1,vjust=0.5),
				strip.placement.x = "inside",
				strip.background=element_blank(),
				panel.spacing.x=unit(0.05, "lines"),
				panel.spacing.y=unit(.2, "lines")
        		)+
		guides(fill=guide_legend(nrow=1,byrow=TRUE, title="Cluster",
									title.position="top"),
				colour=F)



ggsave(filename=out, 
		plot=g, device ="pdf", width=190, height=200, units="mm", dpi ="screen")


system(paste0("open -a Skim.app ", out))
}
plotAdmix()








