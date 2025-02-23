# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-07-13
# Last Modified: 2020-08-14

# Desc: 


###########################################
################# Modules #################
###########################################

require(pophelper)
require(stringr)
require(tools)
require(gridExtra)
require(grid)
require(tidyverse)
require(ggpubr)
require(RColorBrewer)
require(colorspace)

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


admixData <- function(dir, clusters, nAnc){

	# select the files
	files <- list.files(dir, full.names=T)
	files_to_use <- subset(files, file_ext(files) == "qopt")
	files_to_use <- str_sort(files_to_use, numeric = TRUE) # natural sort
	files_to_use <- files_to_use[nAnc]

	# clusters
	pop <- read.delim(clusters, header = TRUE, sep="\t", colClasses = "character")
	# read files to qlist 
	qlist_all <- readQ(files_to_use)
	slist <- alignK(qlist_all, type="auto") # sort "within", "auto", "across"
	k_labs <- c(paste0("K=", nAnc))
	names(slist) <- c(paste0("K", nAnc))

	if(F){

			###
			# manually correct the pops that weren't autocorrected
			# K3 
				# c3 > c2
			names(slist$K3) <- c("Cluster1", "Cluster3", "Cluster2")
			# K4 
				# c1 > c2
				# c4 > c1
			names(slist$K4) <- c("Cluster2", "Cluster4", "Cluster3", "Cluster1")

	}

	# K5
	names(slist$K5) <- c("Cluster2", "Cluster3", "Cluster5", "Cluster4", "Cluster1")
	
	# K6
	names(slist$K6) <- c("Cluster4", "Cluster3", "Cluster5", "Cluster6", "Cluster1", "Cluster2")
	
	# K7
	names(slist$K7) <- c("Cluster1", "Cluster2", "Cluster7", 
								"Cluster4", "Cluster5", "Cluster6", "Cluster3")
	# K8
	names(slist$K8) <- c("Cluster1", "Cluster2", "Cluster7", 
								"Cluster4", "Cluster5", "Cluster6", "Cluster3", "Cluster8")
	# K9
		# c3 ><c7
	names(slist$K9) <- c("Cluster1", "Cluster2", "Cluster7", 
								"Cluster4", "Cluster5", "Cluster6", 
								"Cluster3", "Cluster8", "Cluster9")

	# K15
	names(slist$K15) <- c("Cluster6", "Cluster14", "Cluster2", 
								"Cluster4", "Cluster1", "Cluster10", 
								"Cluster9", "Cluster8", "Cluster15",
								"Cluster5", "Cluster11", "Cluster12",
								"Cluster13", "Cluster3", "Cluster7")

	# K20
	names(slist$K20) <- c("Cluster1", "Cluster3", "Cluster17", 
								"Cluster4", "Cluster16", "Cluster9", 
								"Cluster19", "Cluster13", "Cluster6",
								"Cluster7", "Cluster11", "Cluster12",
								"Cluster8", "Cluster14", "Cluster15",
								"Cluster5", "Cluster20", "Cluster18",
								"Cluster10", "Cluster2")
	# K25
	names(slist$K25) <- c("Cluster5", "Cluster2", "Cluster14", 
								"Cluster10", "Cluster13", "Cluster20", 
								"Cluster4", "Cluster8", "Cluster9",
								"Cluster7", "Cluster11", "Cluster12",
								"Cluster1", "Cluster3", "Cluster15",
								"Cluster16", "Cluster17", "Cluster18",
								"Cluster19", "Cluster6", "Cluster21", 
								"Cluster22", "Cluster23","Cluster24", 
								"Cluster25")

	# K30
	names(slist$K30) <- c("Cluster7", "Cluster2", "Cluster25", 
								"Cluster4", "Cluster29", "Cluster16", 
								"Cluster21", "Cluster8", "Cluster30",
								"Cluster13", "Cluster11", "Cluster3",
								"Cluster10", "Cluster14", "Cluster15",
								"Cluster6", "Cluster17", "Cluster18",
								"Cluster19", "Cluster20", "Cluster1", 
								"Cluster22", "Cluster23","Cluster24", 
								"Cluster5", "Cluster26", "Cluster27", 
								"Cluster12", "Cluster28","Cluster9")

	
	# clusters
	onelabset1 <- pop[,3,drop=FALSE]
	#onelabset1$CLUSTER <- as.vector(sapply(onelabset1$CLUSTER, function(x) substr(x, 1, 3)))
	colnames(onelabset1) <- "grp"
	onelabset1[onelabset1 == "BEN",] <- " BEN "
	onelabset1[onelabset1 == "BENSON",] <- " BENSON "


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

	clstlvls <- as.numeric(unique(store_join$clst))
	store_join$clst <- factor(store_join$clst, levels =clstlvls)

	# create point df to plot
	pointdf <- store_join %>%
				filter(K == "K=2") %>%
				select(grp, val, X, clst) 

	pointdf$val <- 0

	# rename
	store_join$grp <- correctNames(store_join$grp)

	return(store_join)

}


discreteColours <- function(ncol){

	cb <- c("#999999", "#E69F00", "#56B4E9", 
				"#009E73", "#F0E442", "#0072B2", "#D55E00", 
				"#CC79A7", "#000000")
	nleft <- ncol - length(cb)
	colpal <- divergingx_hcl(nleft, palette = "Temps")

	return(c(cb, colpal))

}

plotAdmix <- function(df, out, pwidth=190, pheight=200){


		uniqK <- unique(df$K)
		ncol <- max(sapply(uniqK, function(x) as.numeric(gsub("K=", "", x))))
		largedf <- ncol > 9
		
		if(!largedf){

			colpal <- c("#999999", "#E69F00", "#56B4E9", 
						"#009E73", "#F0E442", "#0072B2", "#D55E00", 
						"#CC79A7", "#000000")
		}else{
			
			#colpal <- wesanderson::wes_palette("Darjeeling1", ncol, type = "continuous")
			colpal <- discreteColours(ncol=ncol)
			
			}


		g <- ggplot(df, aes(fill=clst, colour=clst, y=val, x=X, label=clst))+
				#facet_grid(rows = vars(K))+
				facet_grid(vars(K), vars(grp), 
						scales = "free", space = "free", 
						switch="x") + 
				geom_bar(position=position_stack(reverse = TRUE), 
						stat="identity", width=1,
						size=0.01
						)+
				scale_y_continuous(expand=c(0,0))+
				xlab("")+
				ylab("")+
				scale_fill_manual(values = colpal) +
				scale_colour_manual(values = colpal) +
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
		        		)
				
				if(largedf){
						g <- g + guides(fill=guide_legend(nrow=3,byrow=TRUE, 
									title="Cluster",
									title.position="top"),colour=F)

					}else{
						g <- g + guides(fill=guide_legend(nrow=1,byrow=TRUE, 
									title="Cluster",
									title.position="top"),colour=F)
					}


		ggsave(filename=out, 
				plot=g, device ="pdf", width=pwidth, height=pheight, units="mm", dpi ="screen")


		system(paste0("open -a Skim.app ", out))
}




###########################################
############### plotting #################
###########################################

# all data
dataAll <- admixData(dir="results/ancestry/ALL_5kb_02maf/",  
					clusters="results/ancestry/clusters", 
					nAnc = seq(2, 39))
plotAdmix(df=dataAll, out="results/ancestry/admixture_all.pdf",
					pwidth=190,
					pheight=250)

# selected 
data_select <- admixData(dir="results/ancestry/ALL_5kb_02maf/",  
					clusters="results/ancestry/clusters", 
					nAnc=c(5,6,7,8,9,15,20,25,30))
plotAdmix(df=data_select, out="results/ancestry/admixture.pdf")




