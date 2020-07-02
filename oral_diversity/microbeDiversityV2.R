# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-07-02
# Last Modified: 2020-07-02

# Desc: 

###########################################
################# Modules #################
###########################################

require(ggplot2)
require(dplyr)
require(tidyr)
require(stringr)
require(wesanderson)
require(RColorBrewer)
require(vegan)

# pavian for shiny app
#pavian::runApp(port=5000)

###########################################
############## Function(s) ################
###########################################

readStats <- function(passIn, failIn, lineageFile){

	### read pass stats in ###
	statIn = read.csv(passIn, header = F)
	colnames(statIn) <- c("taxid", "count")
	statIn$count <- as.numeric(statIn$count)

	### read fails in ###
	# sum all fails -  unclassified
	failed <- read.csv(failIn, header = F)
	unclassified <- c(0, sum(failed$V2))

	### read lineages in ###
	lineages <- read.delim(lineageFile, sep="\t", 
						header=T, na.strings = c("","NA"))
	# determine ranks from lineages
	rank <- apply(lineages, 1, 
					function(x) colnames(lineages)[length(x[!is.na(x)])])
	rank <- replace(rank, rank=="taxid", "root")
	lineages_df <- cbind(lineages, rank)

	lineages_df$rank <- as.character(lineages_df$rank)

	df <- statIn %>% 
			rbind(unclassified) %>%
			left_join(lineages_df, by="taxid") %>%
			arrange(rank, desc(count))

	df[df$taxid == 0,]$rank <- "unclassified"

	return(df)

}

###########################################
######### Input(s) and Parameters #########
###########################################


lineageFile = "results/oral_diversity/taxa.df"
passIn = "results/oral_diversity/stats/stat.pass.csv"
failIn = "results/oral_diversity/stats/stat.fail.csv"

###########################################
############### Wraggling #################
###########################################


df <- readStats(passIn, failIn, lineageFile)
tibble(lineages)

df_BacVir <- df[(df$kingdom %in% c("Bacteria", "Viruses"),]

grp_by <- colnames(df)[-c(1,2)]

df_BacVir <- df %>%
				filter(kingdom %in% c("Bacteria", "Viruses") 
						& rank =="species") %>%
				group_by_at(grp_by) %>%
				summarise(count=sum(count), # correct duplicate taxaid - NCBI issues
							taxid = min(taxid)) %>%
				ungroup()


# abundance
tot_rank_frag <- df_BacVir %>% 
					group_by(rank) %>%
					summarise(rank_count = sum(count))

df_abund <- df_BacVir %>%
				left_join(tot_rank_frag, by='rank') %>% 
				mutate(rel_abundance = count/rank_count)

### species only ###

# bacteria
df_species_bact <- df_abund %>%
					filter(rank %in% c("species")) %>%
					filter(kingdom == "Bacteria")


# virus				
df_species_virus <- df_abund %>%
					filter(rank %in% c("species")) %>%
					filter(kingdom == "Viruses")


###########################################
############### Analysis ##################
###########################################

### diversity and effective species number (ESN) ###
# bacteria 
bacteria_counts <- df_species_bact$count
shannon_bact <- diversity(bacteria_counts, index = "shannon")
eff_species_bact <- exp(shannon_bact) # effective number of species

# virus 
virus_counts <- df_species_virus$count
shannon_virus <- diversity(virus_counts, index = "shannon")
eff_species_virus <- exp(shannon_virus) # effective number of species


print("Bacteria...")
print(paste("shannon: ", shannon_bact))
print(paste("ENS: ", eff_species_bact))
print("Virus...")
print(paste("shannon: ", shannon_virus))
print(paste("ENS: ", eff_species_virus))


###########################################
############### Plotting ##################
###########################################

### take must abundent species ###
	# base on the ENS

### bacteria ###
data_head_bact <- df_species_bact %>% 
					arrange(desc(rel_abundance)) %>%
					head(ceiling(eff_species_bact))

#scale_fill_manual(values = pal)+
#pal <- wes_palette("Darjeeling1", nrow(data_head_bact), type = "continuous")

g <- ggplot(data=df_species_bact, aes(x=species, y=count, fill =species)) +
		facet_grid(vars(phylum), scales="free", space = "free")+
		geom_col(colour='black', show.legend = FALSE)+
		theme_classic() +
		ylab("total read match") +
		xlab("") +
		coord_flip()+
		theme(strip.text.y = element_text(size = 15, angle=0),
				text = element_text(size=20))
#options(scipen=10000)
pdf("results/oral_diversity/oralDiv_bacteria.pdf", 15, 20)
print(g)
invisible(dev.off())


### virus ###
data_head_virus <- df_species_virus %>% 
					arrange(desc(rel_abundance)) %>%
					head(ceiling(eff_species_virus))

pal <- wes_palette("Darjeeling1", nrow(data_head_virus), type = "continuous")

p <- ggplot(data=data_head_virus, aes(x=species, y=count, fill =species)) +
		facet_grid(vars(phylum), scales="free", space = "free")+
		geom_col(colour='black', show.legend = FALSE)+
		scale_fill_manual(values = pal)+
		theme_classic() +
		ylab("total read match") +
		xlab("") +
		coord_flip()+
		theme(strip.text.y = element_text(size = 15, angle=0),
				text = element_text(size=20))

pdf("results/oral_diversity/oralDiv_virus.pdf", 15, 20)
print(p)
invisible(dev.off())

system("open -a Skim.app results/oral_diversity/oralDiv_bacteria.pdf")
system("open -a Skim.app results/oral_diversity/oralDiv_virus.pdf")


#########################
### Stats for writeup ###

# number of 
	# plants, archea, fungi, virus

#unique(df_total$kingdom)


#df_total %>%
#	mutate(domain_supking = 
#			ifelse(is.na(kingdom), name, as.character(kingdom))) %>%
#	dplyr::select(-kingdom) %>%
#	dplyr::select(name, domain_supking, phylum, class) %>%
#	unique() %>%
#	group_by(domain_supking) %>%
#	summarise(count=n()) #%>%
#	#write.table("sandbox/fungi_plant.txt", quote=FALSE, row.names = F)


