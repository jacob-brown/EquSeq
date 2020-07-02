# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-15
# Last Modified: 2020-07-01

# Desc: plot treemix results


args = commandArgs(trailingOnly=TRUE)

library(dplyr)
library(ggplot2)

################################
### plot trees

source("../dependancies/treemix-1.13/src/plotting_funcs.R")

singleMig <- function(){
	### single migrations ###
	pdf("tree.pdf", 30, 15)
	par(mfrow=c(1,2))
	plot_tree("tree.out.no", o = NA, cex = 1, disp = 0.003, plus = 0.01, flip = vector(), arrow = 0.05, scale = T, ybar = 0.1, mbar = T, plotmig = T, plotnames = T, xmin = 0, lwd = 1, font = 1)
	plot_resid("tree.out.no", "poporder")
	dev.off()
	system("open -a Skim.app tree.pdf")
}
singleMig()
############################
### multiple  migrations ###
muliMig <- function(){
	prefix="tree.out"
	pdf("tree.pdf", 10, 8)
	par(mfrow=c(2,3))
	for(edge in 0:5){
	  plot_tree(cex=0.8,paste0(prefix,".",edge))
	  title(paste(edge,"edges"))
	}
	dev.off()
	system("open -a Skim.app tree.pdf")

	pdf("residuals.pdf", 20, 15)
	par(mfrow=c(2,3), pin=c(4, 4))
	for(edge in 0:5){
	 plot_resid(stem=paste0(prefix,".",edge),pop_order="poporder")
	}
	dev.off()
	system("open -a Skim.app residuals.pdf")
}
###################################
### check migration likelihoods ###
migLike <- function(){

	files <- list.files("out_treemix", pattern ="*.llik", full.names=T)
	df <- data.frame()
	for(num in 1:length(files)){
		if(length(readLines(file(files[num], "r"))) > 0){ # remove problem repeats
			liketab <- read.table(files[num], sep = ":")
			df <- rbind(df, data.frame((num-1), liketab))
		}
	}
	colnames(df) <- c("mrequest", "output", "likelihood")
	# use final migration state
	df$output <- as.character(df$output)
	con <- apply(df, 2, function(x) grepl("Exiting", x, fixed = TRUE))[,2]
	df_exit <- df[con,]
	df_exit$migrations <- as.factor(gsub("\\D", "", df_exit$output))

	### anova ###
	head(df_exit)
	summary(df_exit)
	mod <- aov(likelihood~migrations, df_exit)
	summary(mod)
	TukeyHSD(mod)

	# for repeats
	pdf("choose_m.pdf", 10, 5)
	par(mfrow=c(1,2))
	boxplot(likelihood~migrations, df_exit)
	plot(TukeyHSD(mod), cex=0.1, las=1)
	dev.off()
	system("open -a Skim.app choose_m.pdf")

}

################################
### f3 stat interpretation
f3 <- function(){

	data <- read.delim("f3stat.txt", stringsAsFactors = FALSE)

	res <- list()
	misc <- list()
	for(i in 1:nrow(data)){
		if(grepl(";", data[i,], fixed=T)){
			res <- append(res, data[i,])
		}else{
			misc <- append(misc, data[i,])
		}
	}

	#head(unlist(res))
	#str = "BENSON;NativeMongolianChakouyiHorse,SwissWarmblood 0.00596735 0.00501543 1.1898" 
	#unlist(strsplit(str, "\\;|\\,| "))

	ls_res_split <- lapply(unlist(res), function(str) unlist(strsplit(str, "\\;|\\,| ")))
	res_split <- data.frame(matrix(unlist(ls_res_split), nrow=length(ls_res_split), byrow=T))
	#head(res_split)
	colnames(res_split) <- c("a", "b", "c", "f3", "stderr", "zscore")


	# modified from: http://popgen.dk/popgen19/pass/slides/Session7_F-stats%20tutorial%20worksheet%20Thursday%20morning.html
	p <- res_split %>%
			filter(b=="BENSON") %>%
	    	mutate(is_significant=zscore <-3) %>%
	    	ggplot(aes(x=a, y=f3, ymin=f3-3*stderr, ymax = f3+3*stderr, color=is_significant)) +
	    	geom_point() + geom_errorbar() + geom_hline(yintercept = 0)

	pdf("f3stat.pdf", 15, 15)
	print(p)
	dev.off()
	system("open -a Skim.app f3stat.pdf")


}

############################

if(args[1] == "s"){
	singleMig() # single population 
}else if(args[1] == "m"){
	muliMig() # multiple migrations 
}else if(args[1] == "l"){
	migLike() # migration liklihood	
}else if(args[1] == "f"){
	f3() # f3 stats
}




########
# ape plotting testing
# cp tree.out.no.treeout.gz tree.ape.gz
# gunzip tree.ape.gz
if(F){

require(ape)

tree <- read.tree(gzfile("tree.out.no.treeout.gz"))
#tree$edge <- read.table(gzfile("tree.out.no.edges.gz"), as.is  = T, comment.char = "", quote = "")


pdf("ape.pdf", 15, 15)
plot(tree)
dev.off()
system("open -a Skim.app ape.pdf")

install.packages("phyloseq")

require(popcorn)
source("../scripts/popcorn_treemix.R")
tree <- read_treemix("tree.out.no")
pdf("pop.pdf", 15, 15)
p <- plot_treemix(tree) + theme_treemix()
print(p)
dev.off()
system("open -a Skim.app pop.pdf")




}