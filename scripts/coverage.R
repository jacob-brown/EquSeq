# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-08-12
# Last Modified: 2020-08-12

# Desc: calculate coverage from samtools coverage results 


###########################################
################# Modules #################
###########################################



###########################################
############## Function(s) ################
###########################################

coverage <- function(fileIn){

	#return average coverage (redundent) and % coverage

	df <- read.delim(fileIn)
	df$X.rname <- as.character(df$X.rname)
	grepl("chr", "crass", fixed = TRUE)
	condition <- as.vector(sapply(df$X.rname, function(x) grepl("chr", x, fixed = TRUE)))
	df <- df[condition, ]
	df <- df[df$X.rname != "chrX" & df$X.rname != "chrM", ]

	# mean coverage (redundent)
	avgcov <- mean(df$meandepth)

	# percentage coverage
	perccov <- (sum(df$covbases) / sum((df$endpos - df$startpos))) * 100

	return(cbind(avgcov, perccov))
}

###########################################
######### Input(s) and Parameters #########
###########################################

all_files <- list.files("results/coverage", full.names =T)

# ensure to only use files that aren't empty
	# error when generating the stats
info <- file.info(all_files)
files <- rownames(info[info$size != 0, ]) 

###########################################
############### Wraggling #################
###########################################

res <- t(sapply(files, function(x) coverage(x)))
df <- data.frame(res)
colnames(df) <- c("avgcov", "perccov")

df_benson <- df[rownames(df)=="results/coverage/coverage_benson.csv", ]
df_all <- df[rownames(df)!="results/coverage/coverage_benson.csv", ]

# return results
cat("Benson")
df_benson
cat("average coverage")
mean(df_all$avgcov)
cat("average % coverage")
mean(df_all$perccov)


###########################################
############### Analysis ##################
###########################################

coverage(fileIn="results/coverage/ERR1545180.sorted.bam.csv")







