# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-07-15
# Last Modified: 2020-08-18

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

correctNames <- function(popls, space=T){
	
	popstore <- as.character(popls)

	correct_df <- read.csv("data/ancestry/clusters_alt.csv", header=F, stringsAsFactors=F)
	colnames(correct_df) <- c("Pop", "new")

	to_change <- correct_df#[correct_df$Pop != correct_df$new,]

	# remove space?
	if(!space){
		to_change$Pop <- as.vector(sapply(to_change$Pop, function(x) gsub(" ", "", x)))
	}

	runs <- length(popstore)
	for(i in 1:runs){

		if(popstore[i] %in% to_change$Pop){
			popstore[i] <- to_change[to_change$Pop == popstore[i],2] # change the name
		}
	}

	return(popstore)
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

bensonf3df <- tibble(f3_df) %>%
	filter(zscore < -2 & A == "BENSON") %>%
	arrange(zscore) %>%
	mutate(f3 = round(f3, 3), stderr = round(stderr, 3), zscore = round(zscore, 3))

bensonf3df$B <- correctNames(bensonf3df$B, space=F)
bensonf3df$C <- correctNames(bensonf3df$C, space=F)

write.csv(bensonf3df, "results/ancestry/f3.benson.csv", row.names = F)

# how many times do the breeds feature?
f3countdf <- data.frame(table(c(bensonf3df$B, bensonf3df$C))) %>% 
				arrange(desc(Freq))
f3countdf

# select 
	# A as outgroup
	# B as Benson (the target)
	# C and D as the pops we want to test
# Przewalski showed signs of mixture

# use f3 to get outgroup
tibble(f3_df) %>%
				filter(A == "BENSON") %>%
				arrange(desc(zscore)) %>%
				head() %>%
				mutate(f3 = round(f3, 3), 
						stderr = round(stderr, 3), 
						zscore = round(zscore, 3)) %>%
write.csv("results/ancestry/f3.outgroups.benson.csv", row.names = F)



# Przewalski seems best or Przewalski-hybrid

tibble(f4_df) %>%
	filter(zscore < -3 | zscore > 3 ) %>%
	filter(A == "Przewalski" & B == "BENSON" ) %>%
	arrange(zscore)  %>%
	mutate(f4 = round(f4, 3), stderr = round(stderr, 3), zscore = round(zscore, 3)) %>%
	write.csv("results/ancestry/f4_prz.benson.csv", row.names = F)


bensonf4df <- tibble(f4_df) %>%
	filter(zscore < -2 | zscore > 2 ) %>%
	filter(A == "Przewalski" & B == "BENSON" 
		& C != "Przewalski-hybrid"& D != "Przewalski-hybrid"
		) %>%
	arrange(zscore)  %>%
	mutate(f4 = round(f4, 3), stderr = round(stderr, 3), zscore = round(zscore, 3))


# rename
bensonf4df$C <- correctNames(bensonf4df$C, space=F)
bensonf4df$D <- correctNames(bensonf4df$D, space=F)

write.csv(bensonf4df, "results/ancestry/f4_nohybrid.benson.csv", row.names = F)


# Split the table and transform it into a more readable format
	# negative f4  gene flow between B and C
	# positive f4  gene flow between B and D
transdf4df <- bensonf4df %>%
				transmute(f4structure = paste0("(", A,".", B, ";",C, ".",D, ")"), 
							Target = B, 
							mixpop = ifelse(zscore < -2, C,  D), 
							f4 = f4,
							stderr = stderr,
							zscore = zscore
							)
negdf4df <- transdf4df[transdf4df$zscore < -2, ] 
posdf4df <- transdf4df[transdf4df$zscore > 2, ] %>%
				arrange(desc(zscore))


write.csv(negdf4df, "results/ancestry/f4_neg.csv", row.names = F)
write.csv(posdf4df, "results/ancestry/f4_pos.csv", row.names = F)

### combine f4 tables for report
f4_all <- rbind(posdf4df, negdf4df)
f4_all <- f4_all[rev(order(abs(f4_all$zscore))),] # order absolute values
f4_all <- f4_all %>% 
	mutate(f4structure = str_replace_all(f4structure, "Przewalski", "A"),
		f4structure = str_replace_all(f4structure, "BENSON", "B"),
		f4structure = str_replace_all(f4structure, ";", "; ")
		) %>%
			transmute("f4 Structure" = str_replace_all(f4structure, "[.]", ","), 
				"Target Mixed With" = mixpop,
				"f4 value" = f4,
				SE=stderr,
				"Z-Score" = zscore)

write.csv(f4_all, "results/ancestry/f4_benson.csv", row.names = F)

# how many times do the breeds feature?
f4countdf <- data.frame(table(c(transdf4df$mixpop))) %>% 
				arrange(desc(Freq))
f4countdf



	