# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-15
# Last Modified: 2020-07-10

# Desc: cluster file from vcf for treemix


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
  	vret <- vapply(strsplit(toupper(charVec), " "), 
                     function(x) paste(substr(x, 1, 1), collapse = ""), 
                     	vector("character", 1L))#)
  	vret[vret == "NA"] <- NA # correct NAs
  	return(vret)
}

###########################################
######### Input(s) and Parameters #########
###########################################
vcffile <- as.character(args[1])
#vcffile <- as.character("snps.rename.vcf")
runName <- system(sprintf("bcftools query -l %s", vcffile), intern = T)
info <- read.csv(args[2])
#info <- read.csv("../data/cleaned_data/info_all.csv")
#info <- read.csv("/rds/general/user/jb1919/home/genomics/EquSeq/data/cleaned_data/info_all.csv")
out <- args[3]
# out <- "ALL.clst"

###########################################
############### Wraggling #################
###########################################

#runName <- apply(files, 1, base)
runDf <- data.frame(index = seq(1,length(runName)), name = runName)

# match the code with info_all
info_trim <- info[c("Run", "sub_group", "BioSample")]

# join on runs
run_join <- left_join(runDf, info_trim, by=c("name" = "Run"))
#head(run_join)
# join on biosamples
biosam_df <- run_join[is.na(run_join$sub_group),1:2]
biosam_names <- unique(biosam_df$name)

# unique df of biosample and cluster name
biosam_clust <- unique(info_trim[info_trim$BioSample %in% biosam_names,2:3])
run_join_biosam <- left_join(biosam_df, biosam_clust, by=c("name" = "BioSample"))

# remove nas from run_join and rbind run_join_biosam
df_all <- run_join %>%
				filter(!is.na(sub_group)) %>%
				dplyr::select(-BioSample) %>%
				rbind(run_join_biosam) %>%
				arrange(index)


#####
# further improve cluster names
# make appreviations
#df_all$sub_group <- makeInitials(df_all$sub_group)

df_all$sub_group <- as.character(df_all$sub_group)

# add benson
df_all[is.na(df_all$sub_group),]$sub_group <- "BENSON"

# correct some of the cluster names
df_all$sub_group <- gsub("\\(>50%Quarter)","50Quarter", gsub(" ", "", df_all$sub_group))


# create and write table
df_all$name <- as.character(df_all$name)

table <- cbind(df_all$name, df_all$name, df_all$sub_group)
#table <- cbind(df_all$name,df_all$name, apprev)
write.table(table, row.names=F, sep="\t", file=out, quote=F, col.names=F)



#clst  <- c(rep("pop1", 10), rep("pop2", 10), rep("pop3", 10), rep("pop4", 15))
#table <- cbind(df_all$name,df_all$name, clst)