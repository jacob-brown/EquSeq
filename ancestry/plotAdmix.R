# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-05-11
# Last Modified: 2020-05-29

# Desc: 

# Rscript scripts/plotAdmix.R -i results/ancestry/ALL.MIX.qopt -o  results/ancestry/ALL.MIX.pdf &> /dev/null

###########################################
################# Modules #################
###########################################

require(ggplot2)
require(tools)
require(dplyr)
require(RColorBrewer)

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


admixture <- function(fileIn, orderIndex){

	file_read_in <- read.table(fileIn)

	### sort by order, then row ###
		# corrects random order of population assignment
	file_ls <- file_read_in[order(match(rownames(file_read_in), orderIndex)),]

	# row - source as reference 
	#ref_trans <- t(file_ls[10,])
	####### CHANGE HERE file_ls[10,] ########
	#ref_df <- data.frame(ref_trans, rownames(ref_trans))
	#colnames(ref_df) <- c("value", "pop")
	#ref_order <- arrange(ref_df, desc(value))
	#col_order <- as.vector(ref_order$pop)
	#file <- file_ls[, col_order]
	file <- file_ls
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

DIR <- "results/ancestry/wg_5kb_02maf/"

files <- list.files(DIR, full.names=T)
out <- "results/ancestry/ALL.MIX.pdf"
clusters <- "results/ancestry/clusters"

cbPalette <- c( "#999999", "#E69F00", "#56B4E9", 
				"#009E73", "#F0E442", "#0072B2", "#D55E00", 
				"#CC79A7", "#000000")

###########################################
############### Wraggling #################
###########################################

# select the files
files_to_use <- subset(files, file_ext(files) == "qopt")

# one log file for sites used
logs <- subset(files, file_ext(files) == "log")[1]
d <- read.table(logs, sep="\n")
sites <- strsplit(as.character(d[6,]), " ")[[1]][7]
maf <- strsplit(as.character(d[4,]), " ")[[1]][3]
title <- paste(sites, " ", maf)

files_to_use <- files_to_use[files_to_use != paste0(DIR,"/ALL.MIX.K10.qopt") &
files_to_use != paste0(DIR,"/ALL.MIX.K11.qopt")][1:5]


i <- 0
store_admix <- 0
orderIndex <- idOrder(files_to_use[2])
orderIndex <- seq(1,45)
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

# plots differ in the groups assigned

#cbPalette <- brewer.pal(K, "Paired")
#cbPalette <- brewer.pal(K, "BrBG")
#admixture(files_to_use[6], orderIndex) %>% filter(CLUSTER == "BENSON")

#cbPalette[1:5]
#cbPalette <- sample(cbPalette, length(cbPalette))

# rename
names <- colnames(store_admix)
names[names=="ind"] <- "pop"
colnames(store_admix) <- names

store_admix$K <- paste0("K=", store_admix$K)

### plot ####
g <- ggplot(store_admix, aes(fill=pop, y=values, x=IID, label=pop))+
		facet_grid(rows = vars(store_admix$K))+
		geom_bar(position=position_stack(reverse = TRUE), stat="identity")+
		scale_fill_manual(values = cbPalette) +
		scale_x_discrete(labels=breed_order)+
		theme_classic()+
		xlab("")+
		ylab("")+
		theme(axis.text.x = element_text(angle = 90, hjust = 1),
				legend.position = "none",
				axis.title.y=element_blank(),
        		axis.text.y=element_blank(),
        		axis.ticks.y=element_blank(),
				strip.text.y = element_text(size = 15, angle=0)
        		)
		#ggtitle(title) +
		#labs(fill = "Ancestral\npopulation")+ 
		#scale_x_discrete(labels=breed_order)+

pdf(file=out, 6, 6)
print(g)
invisible(dev.off())




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