# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-07-15
# Last Modified: 2020-07-15

# Desc: 


###########################################
################# Modules #################
###########################################

require(tidyverse)

###########################################
############## Function(s) ################
###########################################


fstatToDF <- function(fileIn){

	cat(paste0("reading: ", fileIn))

	data <- read.delim(fileIn, stringsAsFactors = FALSE, header=F)

	# remove unwanted headers
	blacklist <- "npop|Estimating|total_nsnp"
	
	tmp <- data %>%
				filter(!grepl(blacklist, V1)) %>%
				mutate(V1 = str_replace_all(V1, ";|\\s|,", "|")) 
	
	npops <- length(strsplit(tmp[1,1], "\\|")[[1]]) - 3

	header <- c(LETTERS[seq(1:npops)], "f", "stderr", "zscore")

	res_split <- separate(tmp, col = V1, sep = "\\|", into = header)

	# convert format
	res_split$f <- as.numeric(res_split$f)
	res_split$stderr <- as.numeric(res_split$stderr) 
	res_split$zscore <- as.numeric(res_split$zscore) 
	
	# rename f3/4
	header[header=="f"] <- paste0("f",npops)
	colnames(res_split) <- header

	return(res_split)

}


###########################################
######### Input(s) and Parameters #########
###########################################

f3File = "results/ancestry/fstat/f3stat.txt"
f4File = "results/ancestry/fstat/f4stat.txt"
fileOut = "results/ancestry/f3.all.csv"
###########################################
############### Wraggling #################
###########################################

f3_df <- fstatToDF(f3File)
f4_df <- fstatToDF(f4File)

tibble(f3_df) %>%
	filter(zscore < -2 & A == "BENSON") %>%
	arrange(zscore) %>%
	write.csv("results/ancestry/f3.benson.csv", row.names = F)

tibble(f4_df) %>%
	filter(zscore < -3 | zscore > 3 ) %>%
	filter(A == "BENSON") %>%
	arrange(desc(zscore))  %>%
	write.csv("results/ancestry/f4.benson.csv", row.names = F)




###########################################
############### Analysis ##################
###########################################



###########################################
############### Plotting ##################
###########################################

	