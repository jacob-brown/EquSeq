# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-03-18
# Last Modified: 2020-05-19

# Desc: 


###########################################
################# Modules #################
###########################################

require(ggplot2)
require(dplyr)
require(tidyr)
require(wesanderson)

# pavian for shiny app
#pavian::runApp(port=5000)

###########################################
############## Function(s) ################
###########################################

### standard error ###
se <- function(x){
	sd(x)/(sqrt(length(x)))
}

### summary of taxa file in ###
tableRank <- function(fileIn){

	rank <- read.delim(fileIn, sep="\t", header=F)
	un <- read.delim('results/oral_diversity/unclassified.txt', sep="\t", header=F)
	root <- read.delim('results/oral_diversity/root.txt', sep="\t", header=F)
	header <- c("tableRank", "nFrag", "nOnlyFrag", "rank", "ncbiTaxa", "name", "file")
	colnames(rank) <- header

	un_root <- cbind(un[c(1,2)], root[c(1,2,7)])
	colnames(un_root) <- c("unPerc", "unFrag", "roPerc", "roFrag", "file")
	totalFrag <- un_root$unFrag + un_root$roFrag 
	un_root_tot <- cbind(un_root, totalFrag) 

	grp_rank <- rank %>% 
					left_join(un_root_tot, by="file") %>%
					group_by(name) %>%
					summarise(avgCov = mean(nFrag), sd = sd(nFrag), 
								se=se(nFrag), avgFragTot = mean(totalFrag)) %>%
					mutate(perc = (avgCov/avgFragTot)*100) %>%
					arrange(desc(perc))

	return(grp_rank)

}


###########################################
######### Input(s) and Parameters #########
###########################################


domain <- tableRank('results/oral_diversity/domain.txt')
phylum <- tableRank('results/oral_diversity/phylum.txt')
genus <- tableRank('results/oral_diversity/genus.txt')
species <- tableRank('results/oral_diversity/species.txt')

###########################################
############### Wraggling #################
###########################################

df <- head(species, 20)

pal <- wes_palette("Darjeeling1", nrow(df), type = "continuous")

g <- ggplot(data=df, aes(x=name, y=perc, fill =name)) +
		geom_col(colour='black', show.legend = FALSE)+
		scale_fill_manual(values = pal)+
		coord_flip() +
		theme_classic() +
		ylab("Percentage Detected") +
		xlab("") #+
		#ylim(0, 60)

pdf("results/oral_diversity/oralDiv.pdf", 5, 5)
print(g)
invisible(dev.off())
system("open -a Skim.app results/oral_diversity/oralDiv.pdf")



###########################################
############### Analysis ##################
###########################################




###########################################
############### Plotting ##################
###########################################



