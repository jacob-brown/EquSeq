# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-05-11
# Last Modified: 2020-08-18

# Desc: 


###########################################
################# Modules #################
###########################################

require(dplyr)

###########################################
############## Function(s) ################
###########################################



###########################################
######### Input(s) and Parameters #########
###########################################

files <- read.table("data/ancestry/bam.list")
info <- read.csv("data/cleaned_data/info_all.csv")

###########################################
############### Wraggling #################
###########################################


#files = rbind(files, files_ww) # remove me later

# strip the path and extension
base <- function(x) regmatches(basename(x), regexpr("^([^.]+)", basename(x))) 

runName <- apply(files, 1, base)

runDf <- data.frame(index = seq(1,length(runName)), name = runName)

# match the code with info_all
info_trim <- info[c("Run", "BioSample", "sub_group", "era")]

run_join <- left_join(runDf, info_trim, by=c("name" = "Run"))

# if Run doesn't match, possible that merging has occured
	# match to BioProject 

run_join$sub_group <- as.character(run_join$sub_group)
run_join$sub_group[run_join$name == "final"] <- "BENSON"
run_join$era[run_join$name == "final"] <- "modern"


# strip na groups for alternative join on biosample
if(any(is.na(run_join$sub_group))){
	run_join_true <- run_join[which(!is.na(run_join$sub_group)),]
	no_group <- subset(run_join, is.na(run_join$sub_group))[,1:2]
	info_uniq <- unique(info_trim[,2:4]) # group by biosample due to multiple runs
	no_joined <- left_join(no_group, info_uniq, by=c("name" = "BioSample"))
	run_join_true <- run_join_true[-3]
	run_join <- rbind(run_join_true, no_joined)

}

# pop "BioSample" if exists
if("BioSample" %in% colnames(run_join)){
	run_join = run_join[which(colnames(run_join)!="BioSample")]
}

# check if other na's have passed
if(any(is.na(run_join$sub_group))){print("check bam.list and info file, NA found")}

	# then final.bam
	# if still no match, output error
len <- length(run_join$sub_group)

table <- cbind(run_join$index,rep(1,len), run_join$sub_group)
colnames(table) <- c("IID","FID","CLUSTER")
df <- data.frame(table)
df$IID <- as.numeric(as.character(df$IID))
df <- df[order(df$IID),]

# write table out
write.table(df, row.names=F, sep="\t", file="results/ancestry/clusters", quote=F)



######
# Summary tables for report writeup

##################
### same a table of the individuals used and accession codes
df_info_all <- run_join %>% 
				filter(index != 90 & sub_group != "BENSON") %>%
				transmute(ascension_biosample = name,
							breed = sub_group
							) %>%
				arrange(breed, ascension_biosample)

write.csv(df_info_all, "results/individs_used.csv", row.names=F, quote=F)


### find the number of each breed used
breed_grp <- read.csv("data/ancestry/breed_grps.csv")
#rbind(data.frame(breed="BENSON", n=1)) %>%
breed_grp_update <- df_info_all %>%
						group_by(breed) %>%
						summarise(n=n()) %>%
						left_join(breed_grp, by="breed") %>%
						arrange(breed) %>%
						transmute(Breed=breed, "Major Type" = major.grp, n)

write.csv(breed_grp_update, "results/breed_grp_update.csv", row.names=F, quote=F)






