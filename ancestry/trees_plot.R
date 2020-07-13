# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-15
# Last Modified: 2020-07-13

# Desc: plot treemix results


args <- commandArgs(trailingOnly=TRUE)

library(dplyr)
library(ggplot2)
library(stringr)

source("dependancies/treemix-1.13/src/plotting_funcs.R")
#prefix <- "results/ancestry/treemix/no_jackknife/tree.out"
#prefix <- "results/ancestry/treemix/no_jackknife/tree.benson.out"
#dir_jack <- "results/ancestry/treemix/jackknife/"
out <- "results/ancestry/"

#pop_order = "results/ancestry/treemix/poporder.dot"
# poporder should use . not -
#stem ="results/ancestry/treemix/no_jackknife/tree.out.5"


################################
### plot trees
singleMig <- function(out){
	
	prefix <- "results/ancestry/treemix/tree.out.0"
	outtree <- paste0(out, "/tree.pdf")
	poporder <- paste0(out, "/treemix/poporder")

	### single migrations ###
	pdf(outtree, 30, 15)
	options(scipen=5)
	par(mfrow=c(1,2))
	 invisible(capture.output(x <- 
	  	plot_tree(prefix, plus=0.0001, disp=0.00003, cex=0.8)
	  	))
	invisible(capture.output(x <- 
		plot_resid(prefix, poporder)
		))
	dev.off()
	#system(paste0("open -a Skim.app ", outtree))
}

############################
### multiple  migrations ###
muliMig <- function(prefix, out){
	outtree <- paste0(out, "/tree.migr.pdf")
	outres <- paste0(out, "/tree.residuals.pdf")
	poporder <- paste0(out, "/treemix/poporder")

	pdf(outtree, 20, 15)
	par(mfrow=c(2,3))
	options(scipen=5)
	for(edge in 0:5){
	  # suppress function output
	  invisible(capture.output(x <- 
	  	plot_tree(paste0(prefix,".",edge), plus=0.0001, disp=0.00003, cex=0.8)
	  	))
	  title(paste(edge,"edges"))
	}
	dev.off()
	system(paste0("open -a Skim.app ", outtree))

	pdf(outres, 20, 15)
	par(mfrow=c(2,3), pin=c(4, 4))
	for(edge in 0:5){
	 invisible(capture.output(x <- 
	 	plot_resid(stem=paste0(prefix,".",edge),pop_order=poporder)
	 	))
	}
	dev.off()
	system(paste0("open -a Skim.app ", outres))
}
#muliMig(prefix, out)
###################################
### check migration likelihoods ###

migLike <- function(dir, out){
	plotout <- paste0(out, "/tree.choose_m.pdf")

	files <- list.files(dir, pattern ="*.llik", full.names=T)
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

	#head(df_exit)
	summary(df_exit)
	mod <- aov(likelihood~migrations, df_exit)
	summary(mod)
	TukeyHSD(mod)

	# for repeats
	pdf(plotout, 10, 5)
	par(mfrow=c(1,2))
	boxplot(likelihood~migrations, df_exit)
	plot(TukeyHSD(mod), cex=0.1, las=1)
	dev.off()
	system(paste0("open -a Skim.app ", plotout))

}


#migLike(migLike, out)

################################
### Variance per migration
variance <- function(prefix){
	df <- data.frame()
	for(i in 0:5){
		stem <- paste0(prefix, ".", i)
		df <- rbind(df, c(i, get_f(stem)))
	}
	colnames(df) <- c("migration", "variance")
	return(df)
}
################################
### Investigate p values for migration edges
#dir_jack
#if(F){
#files <- list.files(dir_jack, pattern ="*.treeout.gz", full.names=T)
#df_store <- data.frame()
#for(i in 100:length(files)){
#
#mN <- as.factor(str_extract(files[i], pattern = "[0-9]+")) # migration N
#lines <- readLines((gzfile(files[i])))
#if(length(lines) > 2){
#	migonly <- lines[3:length(lines)]
#	print(migonly)
#	}else{
#		print("no mig")
#	}
#}
#
#colnames(df_store) <- c("edge_weight", "jack_weight", "jack_se", "pval", "migN")
#head(df_store)
#}


################################
### f3 stat interpretation

f3 <- function(fileIn, plotout){
	data <- read.delim(fileIn, stringsAsFactors = FALSE)
	res <- list()
	misc <- list()
	for(i in 1:nrow(data)){
		if(grepl(";", data[i,], fixed=T)){
			res <- append(res, data[i,])
		}else{
			misc <- append(misc, data[i,])
		}
	}

	ls_res_split <- lapply(unlist(res), function(str) unlist(strsplit(str, "\\;|\\,| ")))
	res_split <- data.frame(matrix(unlist(ls_res_split), nrow=length(ls_res_split), byrow=T))

	colnames(res_split) <- c("target", "source_1", "source_2", "f3", "stderr", "zscore")

	cat("writing to flat table: results/ancestry/f3.all.csv")
	write.csv(res_split, "results/ancestry/f3.all.csv", row.names = F)

	res_benson <- res_split %>% 
					filter(target=="BENSON") %>%
					mutate(f3 = as.numeric(as.character(f3)))
	write.csv(res_benson, "results/ancestry/f3.benson.csv", row.names = F)

	
	pdf(plotout, 10, 5)
	boxplot(f3~source_1, res_benson)
	dev.off()
	system(paste0("open -a Skim.app ", plotout))

	# modified from: http://popgen.dk/popgen19/pass/slides/Session7_F-stats%20tutorial%20worksheet%20Thursday%20morning.html
	#p <- res_split %>%
	#		filter(target=="BENSON") %>%
	#    	mutate(is_significant=zscore <-3) %>%
	#    	ggplot(aes(x=a, y=f3, ymin=f3-3*stderr, ymax = f3+3*stderr, color=is_significant)) +
	#    	geom_point() + geom_errorbar() + geom_hline(yintercept = 0)
#
#	#pdf(plotout, 15, 15)
#	#print(p)
#	#dev.off()
	#system(paste0("open -a Skim.app ", plotout))


}
#f3(fileIn, plotout)
############################
########## Main ############

if(args[1] == "s"){
	singleMig(out) # single population 
}else if(args[1] == "m"){
	muliMig(prefix, out) # multiple migrations 
}else if(args[1] == "l"){
	migLike(dir_jack, out) # migration liklihood	
}else if(args[1] == "f"){
	f3(fileIn = "results/ancestry/fstat/f3stat.txt", 
		plotout = "results/ancestry/f3.pdf") # f3 stats
}else if(args[1] == "v"){
	variance(prefix)
}


