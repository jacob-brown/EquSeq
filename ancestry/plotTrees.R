# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-15
# Last Modified: 2020-08-17

# Desc: plot treemix results


###########################################
################# Modules #################
###########################################

library(dplyr)
library(ggplot2)

#source("dependancies/treemix-1.13/src/plotting_funcs.R")
source("scripts/jb_plotting_funcs.R")
###########################################
############## Function(s) ################
###########################################

### plot trees
plotTree <- function(outfile, poporder){
	
	prefix <- "results/ancestry/treemix/tree.out.0"

	### single migrations ###
	pdf(outfile, width=6, height=6)
	options(scipen=5)
	 invisible(capture.output(x <- 
	  	plot_tree(stem = prefix, plus=0.0003, disp=0.00001, cex=0.7,
	  		o = poporder, 
	  		flip = vector(), arrow = 0.05, 
	  		scale = T,  # SE
	  		ybar = 0.1, 
	  		mbar = F, plotmig = F, # migration
	  		plotnames = T, xmin = 0, lwd = 1.5, font = 1
	  		)
	  	))
	dev.off()
	system(paste0("open -a Skim.app ", outfile))
}
plotTree(outfile="results/ancestry/tree.pdf",
			poporder = "results/ancestry/treemix/poporder_col")
### plot tree residuals
plotRes <- function(outfile, poporder){
	
	prefix <- "results/ancestry/treemix/tree.out.0"

	### single migrations ###
	pdf(outfile, width=7, height=7)
	par(mar = c(10, 10, 1, 1))
	invisible(capture.output(x <- 
		plot_resid(prefix, poporder,  
			min = -0.009, max = 0.009, 
			cex = 0.7, usemax = T, wcols = "r")

		))
	dev.off()
	#system(paste0("open -a Skim.app ", outfile))
}

strip_under <- function(vec) sapply(vec, function(x) sub("_",".", x))


###########################################
################## Plot ###################
###########################################

# update poporder file - for colouring plot
	# cp poporder  poporder_col
cbPalette <- c( "#000000", "#F0E442", "#E69F00", 
				"#009E73",  
				"#CC79A7",   "#D55E00", "#999999",  "#56B4E9", "#0072B2")

# major breeds
maj_breed <- read.table("data/ancestry/breed_grps_tree.csv", sep = ",", header=T)
maj_breed$breed <- sapply(maj_breed$breed, function(x) gsub(" ", "", x))
maj_breed$major.grp.num <- as.numeric(maj_breed$major.grp)
maj_breed$col <- NA
for(i in unique(maj_breed$major.grp.num)){
	maj_breed[maj_breed$major.grp.num == i, ]$col <- cbPalette[i]
}

poporder_col <- maj_breed[c("breed", "col")]
write.table(poporder_col, "results/ancestry/treemix/poporder_col", 
	row.names = FALSE, col.names = FALSE ,quote = FALSE, sep = " ")



##############
#### Plot ####

# tree with coloured groups 
plotTree(outfile="results/ancestry/tree.pdf",
			poporder = "results/ancestry/treemix/poporder_col")

# residuals - black colour
#plotRes(out="results/ancestry/tree_res.pdf", 
#	poporder = "results/ancestry//treemix/poporder" )

# variance 
get_f("results/ancestry/treemix/tree.out.0")











#############################################
##### Working but unused functions ####



#prefix <- "results/ancestry/treemix/no_jackknife/tree.out"
#prefix <- "results/ancestry/treemix/no_jackknife/tree.benson.out"
#dir_jack <- "results/ancestry/treemix/jackknife/"
#pop_order = "results/ancestry/treemix/poporder.dot"
# poporder should use . not -
#stem ="results/ancestry/treemix/no_jackknife/tree.out.5"

if(FALSE){
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

	############################
	########## Main ############

	if(args[1] == "s"){
		singleMig(out) # single population 
	}else if(args[1] == "m"){
		muliMig(prefix, out) # multiple migrations 
	}else if(args[1] == "l"){
		migLike(dir_jack, out) # migration liklihood	
	#}else if(args[1] == "f"){
	#	f3(fileIn = "results/ancestry/fstat/f3stat.txt", 
	#		plotout = "results/ancestry/f3.pdf") # f3 stats
	}else if(args[1] == "v"){
		variance(prefix)
	}


}

