# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-05-05
# Last Modified: 2020-05-05

# Desc: generate a table of names for running the pcangsd script


###########################################
################# Modules #################
###########################################

require(tools)
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

# strip the path and extension
base <- function(x) regmatches(basename(x), regexpr("^([^.]+)", basename(x))) 

runName <- apply(files, 1, base)

runDf <- data.frame(index = seq(1,length(runName)), name = runName)

# match the code with info_all
info_trim <- info[c("Run", "BioSample", "sub_group", "era")]

run_join <-left_join(runDf, info_trim, by=c("name" = "Run"))

# if Run doesn't match, possible that merging has occured
	# match to BioProject 

#if(any(is.na(run_join$sub_group))){
#
#	run_join = left_join(runDf, info_trim, by=c("name" = "BioSample"))
#
#}

run_join$sub_group <- as.character(run_join$sub_group)
run_join$sub_group[run_join$name == "final"] <- "NOVEL"
run_join$era[run_join$name == "final"] <- "modern"


	# then final.bam
	# if still no match, output error
len <- length(run_join$sub_group)

table <- cbind(run_join$index,rep(1,len), run_join$sub_group)
colnames(table) <- c("FID","IID","CLUSTER")
# write table out

#write.table(table, row.names=F, sep="\t", col.names=c("FID","IID","CLUSTER"), file="results/ancestry/test.clst", quote=F)


###########################################
################## Plot ###################
###########################################

#cov <- as.matrix(read.table("results/ancestry/test.cov"))
#e <- eigen(cov)
#
#
#plot(e$vectors[,1:2], col=ID$CLUSTER)
#legend("topleft",fill=1:4,levels(ID$CLUSTER))

###### run PCA script ########
	# Modified Matteo's plotPCA.R 
		# not enough colour palette options

#com <- "Rscript scripts/plotPCA.R -i results/ancestry/test.cov -c 1-2 -a results/ancestry/test.clst -o results/ancestry/test.pca.pdf"

# exec
#system(com)
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Read input file
covar <- read.table("results/ancestry/test.cov", stringsAsFact=F);

# Read annot file
annot <- read.table("results/ancestry/test.clst", sep="\t", header=T); # note that plink cluster files are usually tab-separated

# Parse components to analyze
comp <- as.numeric(strsplit("1-2", "-", fixed=TRUE)[[1]])

# Eigenvalues
eig <- eigen(covar, symm=TRUE);
eig$val <- eig$val/sum(eig$val);
cat(signif(eig$val, digits=3)*100,"\n");

# Plot
PC <- as.data.frame(eig$vectors)
colnames(PC) <- gsub("V", "PC", colnames(PC))
PC$Pop <- factor(annot$CLUSTER)

title <- paste("PC",comp[1]," (",signif(eig$val[comp[1]], digits=3)*100,"%)"," / PC",comp[2]," (",signif(eig$val[comp[2]], digits=3)*100,"%)",sep="",collapse="")

x_axis = paste("PC",comp[1],sep="")
y_axis = paste("PC",comp[2],sep="")


rainCP = rainbow(length(unique(PC$Pop)))
novel_Index <- match(c("NOVEL"),sort(unique(PC$Pop)))

rainCP[novel_Index] <- 	"#000000"

ggplot() + 
	geom_point(data=PC, aes_string(x=x_axis, y=y_axis, color="Pop")) + 
	ggtitle(title) + 
	scale_colour_manual(values=rainCP)
	# + scale_colour_manual(values=cbPalette)
ggsave("results/ancestry/test.pca.pdf")
unlink("Rplots.pdf", force=TRUE)




