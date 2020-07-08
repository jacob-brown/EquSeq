# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-07-02
# Last Modified: 2020-07-03

# Desc: plot and analyse summarised microbial classification data

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


readKReport <- function(fileIn, lineageFile){

	lineages <- read.delim(lineageFile, sep="\t", 
						header=T, na.strings = c("","NA"))
	tmp_rep <- read.delim(fileIn, sep="\t", head=F)
	colnames(tmp_rep) <- c("perc", "frag", "frag_uniq", "rankAb", "taxid", "name")
	# perc = (frag/sum(tmp_rep$frag_uniq))*100

	# remove whitespace
	tmp_rep$name <- trimws(tmp_rep$name, which="both")

	# determine ranks from lineages
	rank <- apply(lineages, 1, 
					function(x) colnames(lineages)[length(x[!is.na(x)])])
	rank <- replace(rank, rank=="taxid", "root")
	lineages_df <- cbind(lineages, rank)

	tmp_rep$rankAb <- as.character(tmp_rep$rankAb)

	# join report and lineages
	df <- tmp_rep %>%
			left_join(lineages_df, by="taxid") %>%
			mutate(
				rank = 
					ifelse(rankAb == "D", "domain",
						ifelse(rankAb == "U", "unclassified",
						# update rank
						ifelse(grepl("\\d", rankAb),
									paste0(rank,"_", str_extract(rankAb, "[[:digit:]]+")),
									as.character(rank))))#,
				#kingdom = ifelse(name != kingdom, name, kingdom)
					# some kingdoms were defined incorretly
				) %>%
			arrange(rank, desc(perc))

	return(df)

}

###########################################
######### Input(s) and Parameters #########
###########################################

lineageFile = "results/oral_diversity/taxa.df"
passIn = "results/oral_diversity/stat.all.pass.csv"
failIn = "results/oral_diversity/stat.all.fail.csv"
report <- "results/oral_diversity/kraken_reports//V300044309_L2_B5GHORlfyRAAAAAAA-517.kreport"

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", 
				"#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

###########################################
############### Wraggling #################
###########################################


df <- readStats(passIn, failIn, lineageFile)

df_BacVir <- df[df$kingdom %in% c("Bacteria", "Viruses"),]

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
bacteria_total <- length(bacteria_counts)
shannon_bact <- diversity(bacteria_counts, index = "shannon")
eff_species_bact <- exp(shannon_bact) # effective number of species

# virus 
virus_counts <- df_species_virus$count
virus_total <- length(virus_counts)
shannon_virus <- diversity(virus_counts, index = "shannon")
eff_species_virus <- exp(shannon_virus) # effective number of species

cat("\n=====================\n")
cat("\nStatistics\n")
cat("\nBacteria...\n")
cat(paste("shannon: ", shannon_bact, "\n"))
cat(paste("ENS: ", eff_species_bact, "\n"))
cat(paste("Total species detected: ", bacteria_total, "\n"))
cat("\nVirus...\n")
cat(paste("shannon: ", shannon_virus, "\n"))
cat(paste("ENS: ", eff_species_virus, "\n"))
cat(paste("Total species detected: ", bacteria_total, "\n"))


###########################################
############### Plotting ##################
###########################################

### take must abundent species ###
	# base on the ENS

### bacteria ###
data_head_bact <- df_species_bact %>% 
					arrange(desc(rel_abundance)) %>%
					#head(100)
					head(ceiling(eff_species_bact))

#ncolours <- length(unique(data_head_bact$phylum))
#pal <- wes_palette("Darjeeling1", ncolours, type = "continuous")
#scale_fill_manual(values = pal)+


g <- ggplot(data=data_head_bact, aes(x=species, y=count, fill =phylum)) +
		facet_grid(vars(phylum), scales="free", space = "free")+
		geom_col(colour='black', show.legend = FALSE)+
		theme_classic() +
		ylab("total read match") +
		xlab("species") +
		scale_y_continuous(labels = scales::comma)+ 
		coord_flip()+
		scale_fill_manual(values = cbPalette)+
		theme(strip.text.y = element_text(size = 15, angle=0),
				text = element_text(size=20))
#options(scipen=10000)
pdf("results/oral_diversity/oralDiv_bacteria.pdf", 15, 20)
print(g)
invisible(dev.off())


### virus ###
data_head_virus <- df_species_virus %>% 
					arrange(desc(rel_abundance)) %>%
					#head(100)
					head(ceiling(eff_species_virus))

pal <- wes_palette("Darjeeling1", nrow(data_head_virus), type = "continuous")

p <- ggplot(data=data_head_virus, aes(x=species, y=count, fill =phylum)) +
		facet_grid(vars(phylum), scales="free", space = "free")+
		geom_col(colour='black', show.legend = FALSE)+
		theme_classic() +
		ylab("total read match") +
		xlab("species") +
		coord_flip()+
		scale_fill_manual(values = cbPalette)+
		theme(strip.text.y = element_text(size = 15, angle=0),
				text = element_text(size=20))

pdf("results/oral_diversity/oralDiv_virus.pdf", 15, 20)
print(p)
invisible(dev.off())

system("open -a Skim.app results/oral_diversity/oralDiv_bacteria.pdf")
system("open -a Skim.app results/oral_diversity/oralDiv_virus.pdf")


### save datatables ###
write.csv(df_species_bact, "results/oral_diversity/df_bacteria.csv",row.names = F, quote=F)
write.csv(df_species_virus, "results/oral_diversity/df_virus.csv",row.names = F, quote=F)
#########################
### Stats for writeup ###
#########################

data_report  <- readKReport(report, lineageFile)
grps <- setdiff(names(data_report), c("perc", "frag", "frag_uniq"))


# use total frag counts across all the files
df_total <- data_report %>% 
					group_by_at(grps) %>%
					summarise(tot_perc = sum(perc), 
								tot_frag = as.numeric(sum(frag)), 
								tot_frag_uniq = sum(frag_uniq),
								count=n()) %>%
					arrange(desc(tot_perc)) %>%
					ungroup()


kings <- df_total %>%
			filter(!is.na(species)) %>%
			select(kingdom, species) %>% 
			unique() %>%
			group_by(kingdom) %>%
			summarise(count=n()) %>%
			transmute(type = kingdom, count)


fungi <- c("Ascomycota", "Basidiomycota", "Microsporidia")
plants <-c("Chlorophyta", "Rhodophyta", "Streptophyta")

df_total %>%
	filter(kingdom == "Eukaryota" & !is.na(species)) %>%
	mutate(phylum = as.character(phylum),
			type = ifelse(phylum %in% fungi, "fungi",
					ifelse(phylum %in% plants, "plant", phylum))) %>%
	select(type, species) %>%
	unique() %>%
	group_by(type) %>%
	summarise(count=n()) %>%
	union(kings) %>%
	write.table("results/oral_diversity/kraken_db_composition.txt", quote=FALSE, row.names = F)

