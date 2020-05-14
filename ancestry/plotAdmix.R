# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-05-11
# Last Modified: 2020-05-14

# Desc: 

# Rscript scripts/plotAdmix.R -i results/ancestry/ALL.MIX.qopt -o  results/ancestry/ALL.MIX.pdf &> /dev/null

###########################################
################# Modules #################
###########################################

require(ggplot2)
require(tools)
require(dplyr)

###########################################
############## Function(s) ################
###########################################

idOrder <- function(fileIn){
	# individual order for plotting from file
	file <- read.table(fileIn)
	admix <- cbind(file, IID=as.integer(rownames(file)))
	pop <- read.table(clusters, header = TRUE, sep="\t")

	admix_pop <- left_join(admix, pop, by=c("IID"))
	admix_pop <- admix_pop %>% arrange(desc(V1))

	source_row <- subset(admix_pop, admix_pop$CLUSTER=="Source")
	all_rows <- subset(admix_pop, admix_pop$CLUSTER!="Source")
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
	ref_trans <- t(file_ls[10,])
	###### CHANGE HERE file_ls[10,] ########
	ref_df <- data.frame(ref_trans, rownames(ref_trans))
	colnames(ref_df) <- c("value", "pop")
	ref_order <- arrange(ref_df, desc(value))
	col_order <- as.vector(ref_order$pop)
	file <- file_ls[, col_order]

	# transform
	colnames(file) <- letters[1:ncol(file)]
	admix <- utils::stack(file)
	admix <- cbind(admix, IID = rep(1:nrow(file)))

	K <- length(unique(admix$ind)) # nrow(admix)
	if (K > 8) stop("Maximum number of colours for admixture proportions is 8.")

	# annotate with clusters
	pop <- read.table(clusters, header = TRUE, sep="\t")
	admix_pop <- left_join(admix, pop, by=c("IID"))
	admix_pop <- cbind(admix_pop, K=K)
	
	return(admix_pop)


}




###########################################
######### Input(s) and Parameters #########
###########################################

files <- list.files("results/ancestry/", full.names=T)
out <- "results/ancestry/ALL.MIX.pdf"
clusters <- "results/ancestry/clusters"

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
				"#F0E442", "#0072B2", "#D55E00", "#CC79A7")

###########################################
############### Wraggling #################
###########################################

# select the files
files_to_use = c()
for(f in files){
	if(file_ext(f) == "qopt"){
		files_to_use <- append(files_to_use, f)
	}
}

i <- 0
store_admix <- 0

for(i in seq_along(files_to_use)){
	
	# retrieve order on first iter
	if(i == 1){
		orderIndex = idOrder(files_to_use[i])
		store_admix <- admixture(files_to_use[i], orderIndex)

	}else{
		tmp <- admixture(files_to_use[i], orderIndex)
		store_admix <- rbind(store_admix, tmp)
	}
}

# sort
#store_admix <- store_admix[with(store_admix, order(K, ind, values)),]

# levels based on first lowest K
K <- min(store_admix$K)
#sub_admix <- store_admix[which(store_admix$K == K),]
#lev <- rev(sub_admix$IID[1:(length(sub_admix$IID)/K)])
#store_admix$IID <- factor(store_admix$IID, levels = lev)

#store_admix <- store_admix %>% filter(IID==45) 
#store_admix <- store_admix %>% filter(K==4) 
#store_admix <- store_admix %>% filter(K==5) 

### sort ###
# plots differ in the groups assigned


g <- ggplot(store_admix, aes(fill=ind, y=values, x=IID, label=ind))+
		facet_grid(rows = vars(store_admix$K))+
		geom_bar(position=position_stack(reverse = TRUE), stat="identity")+
		scale_fill_manual(values = cbPalette) +
		geom_text()+
		#scale_x_discrete(breaks=pop$IID, expand=c(0,0.5))+
		#scale_y_continuous(expand=c(0,0))+
		theme_classic()

pdf(file=out, 13, 10)
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