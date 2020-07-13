# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-05-11
# Last Modified: 2020-07-13

# Desc: 

# matteos: Rscript scripts/plotAdmix.R -i results/ancestry/ALL.MIX.qopt -o  results/ancestry/ALL.MIX.pdf &> /dev/null

###########################################
################# Modules #################
###########################################

require(ggplot2)
require(tools)
require(dplyr)
require(tibble)
require(RColorBrewer)
require(stringr)

###########################################
############## Function(s) ################
###########################################


admixture <- function(fileIn){

	file_read_in <- read.table(fileIn)
	file <- file_read_in

	# transform
	colnames(file) <- seq(1, ncol(file))
	admix <- utils::stack(file)
	admix <- cbind(admix, IID = rep(1:nrow(file)))

	K <- length(unique(admix$ind)) # nrow(admix)
	#if (K > 8) stop("Maximum number of colours for admixture proportions is 8.")

	# annotate with clusters
	pop <- read.table(clusters, header = TRUE, sep="\t")
	admix_pop <- left_join(admix, pop, by=c("IID"))
	admix_pop <- cbind(admix_pop, K=K)
	
	return(admix_pop)


}


###########################################
######### Input(s) and Parameters #########
###########################################

maxAnc <- 9

DIR <- "results/ancestry/ALL_5kb_02maf/"
files <- list.files(DIR, full.names=T)
out <- "results/ancestry/ALL.MIX.pdf"
clusters <- "results/ancestry/clusters"

cbPalette <- c( "#999999", "#E69F00", "#56B4E9", 
				"#009E73", "#F0E442", "#0072B2", "#D55E00", 
				"#CC79A7", "#000000")

cbPalettelrg <- c("#004949","#009292","#ff6db6","#ffb6db",
 					"#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
 					"#920000","#924900","#db6d00","#24ff24","#ffff6d", "#000000")

###########################################
############### Wraggling #################
###########################################

# select the files
files_to_use <- subset(files, file_ext(files) == "qopt")
files_to_use <- str_sort(files_to_use, numeric = TRUE) # natural sort
files_to_use <- files_to_use[2:maxAnc]

# one log file for sites used
logs <- subset(files, file_ext(files) == "log")[3]
d <- read.table(logs, sep="\n")
sites <- strsplit(as.character(d[6,]), " ")[[1]][7]
maf <- strsplit(as.character(d[4,]), " ")[[1]][3]
title <- paste(sites, " ", maf)

# create one df of all admixture data
ls_df <- lapply(files_to_use, function(x) admixture(x))
store_admix <- do.call(rbind, ls_df)

# levels based on first lowest K
K <- max(store_admix$K)
breeds <- filter(store_admix, K==2)$CLUSTER
breed_order <- sort(as.character(unique(breeds)))
breed_order <- c("BENSON", breed_order[!breed_order== "BENSON"])
store_admix$CLUSTER <- factor(store_admix$CLUSTER, levels = breed_order)

# rename
names <- colnames(store_admix)
names[names=="ind"] <- "pop"
colnames(store_admix) <- names

store_admix$K <- paste0("K=", store_admix$K)
store_admix <- tibble(store_admix)

tmp <- store_admix %>%
			filter(K == "K=2" & pop == 1) %>%
			rowid_to_column(var="rowid") %>%
			data.frame()

lab = c()
br = c()
for(i in 1:nrow(tmp)){
	row <- tmp[i,]
	if(!(row$CLUSTER %in% lab)){
		lab <- append(lab, as.character(row$CLUSTER)) 
		br <- append(br, row$IID) 
	}


}


k_levels <- str_sort(unique(store_admix$K), numeric = TRUE)
store_admix$K <- factor(store_admix$K, levels = k_levels)
clst_lvls <- as.character(unique(store_admix$CLUSTER))
store_admix$IID <- factor(store_admix$IID)
store_admix$clust_iid <- paste0(store_admix$CLUSTER,"_",store_admix$IID)
clst_iid_lvls <- sort(unique(store_admix$clust_iid))
clst_iid_lvls <- c(clst_iid_lvls[!clst_iid_lvls == "BENSON_172"], "BENSON_172")
store_admix$clust_iid <- factor(store_admix$clust_iid, levels = clst_iid_lvls)

### plot ####
#scale_fill_manual(values = cbPalettelrg) +
#scale_x_discrete(labels=breed_order)+
g <- ggplot(store_admix, aes(fill=pop, y=values, x=clust_iid, label=pop))+
		facet_grid(rows = vars(store_admix$K))+
		geom_bar(position=position_stack(reverse = TRUE), stat="identity")+
		theme_classic()+
		xlab("")+
		ylab("")+
		scale_fill_manual(values = cbPalette) +
		theme(axis.text.x = element_text(angle = 90, hjust = 1, size =8),
				 axis.line.y = element_blank(),
				legend.position = "none",
				axis.title.y = element_blank(),
        		axis.text.y = element_blank(),
        		axis.ticks.y = element_blank(),
				strip.text.y = element_text(size = 4, angle=0),
				strip.background=element_blank()
        		)
pdf(file=out, 20, 10)
print(g)
invisible(dev.off())
system("open -a Skim.app results/ancestry/ALL.MIX.pdf")


###########################################
################## Plot ###################
###########################################
# Modified Matteo's plotPCA.R 
	# not enough colour palette options
	# and fine tune ggplot