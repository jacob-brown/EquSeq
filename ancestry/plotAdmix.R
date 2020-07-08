# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-05-11
# Last Modified: 2020-07-08

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

idOrder <- function(fileIn){
	# individual order for plotting from file
	file <- read.table(fileIn)
	admix <- cbind(file, IID=as.integer(rownames(file)))
	pop <- read.table(clusters, header = TRUE, sep="\t")

	admix_pop <- left_join(admix, pop, by=c("IID"))
	# get the highest value for source and use it to order
	source_index <- which(admix_pop$CLUSTER == "BENSON")
	vals <- admix_pop[source_index,1:2]
	
	if(vals[1] > vals[2]){
			admix_pop <- admix_pop %>% arrange(desc(V1))
		}else{
			admix_pop <- admix_pop %>% arrange(desc(V2))
		}


	#admix_pop <- admix_pop %>% arrange(desc(V1))

	source_row <- subset(admix_pop, admix_pop$CLUSTER=="BENSON")
	all_rows <- subset(admix_pop, admix_pop$CLUSTER!="BENSON")
	admix_ordered <- rbind(source_row, all_rows)
	order_iid <- admix_ordered$IID
	return(order_iid)
}


admixtureV1 <- function(fileIn, orderIndex){

	file_read_in <- read.table(fileIn)

	### sort by order, then row ###
		# corrects random order of population assignment
	file <- file_read_in[order(match(rownames(file_read_in), orderIndex)),]
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

admixture <- function(fileIn, orderIndex){

	file_read_in <- read.table(fileIn)
	file <- file_read_in
	### sort by order, then row ###
		# corrects random order of population assignment
	#file <- file_read_in[order(match(rownames(file_read_in), orderIndex)),]
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

maxAnc <- 6

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
files_to_use <- files_to_use[1:maxAnc-1]


# one log file for sites used
logs <- subset(files, file_ext(files) == "log")[3]
d <- read.table(logs, sep="\n")
sites <- strsplit(as.character(d[6,]), " ")[[1]][7]
maf <- strsplit(as.character(d[4,]), " ")[[1]][3]
title <- paste(sites, " ", maf)


store_admix <- 0
orderIndex <- idOrder(files_to_use[2])
#orderIndex <- seq(1,length(orderIndex))
for(i in seq_along(files_to_use)){
	
	# retrieve order on first iter
	if(i == 1){
		#orderIndex = idOrder(files_to_use[i])
		# assign order to global for use in plotting
		store_admix <- admixture(files_to_use[i], orderIndex)


	}else{
		tmp <- admixture(files_to_use[i], orderIndex)
		store_admix <- rbind(store_admix, tmp)
	}
}



# levels based on first lowest K
K <- max(store_admix$K)
breeds <- filter(store_admix, K==2)$CLUSTER
breed_order <- breeds[orderIndex]
store_admix$IID <- factor(store_admix$IID, levels=orderIndex)


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
#store_admix$CLUSTER <- factor(store_admix$CLUSTER, levels = clst_lvls)


### plot ####
g <- ggplot(store_admix, aes(fill=pop, y=values, x=IID, label=pop))+
		facet_grid(rows = vars(store_admix$K))+
		geom_bar(position=position_stack(reverse = TRUE), stat="identity")+
		#scale_fill_manual(values = cbPalettelrg) +
		scale_x_discrete(labels=breed_order)+
		theme_classic()+
		xlab("")+
		ylab("")+
		theme(axis.text.x = element_text(angle = 90, hjust = 1, size =8),
				 axis.line.y = element_blank(),
				legend.position = "none",
				axis.title.y = element_blank(),
        		axis.text.y = element_blank(),
        		axis.ticks.y = element_blank(),
				strip.text.y = element_text(size = 4, angle=0),
				strip.background=element_blank()
        		)
pdf(file=out, 20, 20)
print(g)
invisible(dev.off())


#scale_x_discrete(breaks=br, labels=lab)+
#ggtitle(title) +
#labs(fill = "Ancestral\npopulation")+ 
#scale_x_discrete(labels=breed_order)+


#pop<-read.table(clusters, header = TRUE, sep="\t")
##admix<-t(as.matrix(read.table("myoutfiles.qopt")))
#admix<-admix[,order(pop[,1])]
#pop<-pop[order(pop[,1]),]


#pdf(file=out, 13, 5)
#h<-barplot(admix,col=1:3,space=0,border=NA,xlab="Individuals",ylab="admixture")
#text(tapply(1:nrow(pop),pop[,1],mean),-0.05,unique(pop[,1]),xpd=T)
#
#invisible(dev.off())




###########################################
################## Plot ###################
###########################################
# Modified Matteo's plotPCA.R 
	# not enough colour palette options
	# and fine tune ggplot