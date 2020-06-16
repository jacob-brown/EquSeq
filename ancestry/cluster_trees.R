# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-15
# Last Modified: 2020-06-16

# Desc: 


###########################################
################# Modules #################
###########################################

library(dplyr)

###########################################
############## Function(s) ################
###########################################


args = commandArgs(trailingOnly=TRUE)

makeInitials <- function(charVec) {
  #make.unique(
  	vapply(strsplit(toupper(charVec), " "), 
                     function(x) paste(substr(x, 1, 1), collapse = ""), 
                     vector("character", 1L))#)
}

###########################################
######### Input(s) and Parameters #########
###########################################
vcffile <- as.character(args[1])
runName <- system(sprintf("bcftools query -l %s", vcffile), intern = T)
info <- read.csv("../data/cleaned_data/info_all.csv")

###########################################
############### Wraggling #################
###########################################


#runName <- apply(files, 1, base)
runDf <- data.frame(index = seq(1,length(runName)), name = runName)

# match the code with info_all
info_trim <- info[c("Run", "sub_group")]

run_join <- left_join(runDf, info_trim, by=c("name" = "Run"))

# if Run doesn't match, possible that merging has occured
	# match to BioProject 

run_join$sub_group <- as.character(run_join$sub_group)
run_join$sub_group <- gsub("\\(>50%Quarter)","", gsub(" ", "", run_join$sub_group))
#run_join$sub_group <- gsub(" \\(>50% Quarter)","", run_join$sub_group)

# add benson
run_join[is.na(run_join$sub_group),]$sub_group <- "BENSON"

# make appreviations
#apprev <- makeInitials(run_join$sub_group)

# create and write table



table <- cbind(run_join$name,run_join$name, run_join$sub_group)
#table <- cbind(run_join$name,run_join$name, apprev)
write.table(table, row.names=F, sep="\t", file="clusters.clst", quote=F, col.names=F)



#clst  <- c(rep("pop1", 10), rep("pop2", 10), rep("pop3", 10), rep("pop4", 15))
#table <- cbind(run_join$name,run_join$name, clst)