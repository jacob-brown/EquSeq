# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-03-18
# Last Modified: 2020-06-12

# Desc: Read kraken reports and generate plots and 
	# calculate diversity


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

### standard error ###
#se <- function(x){
#	sd(x)/(sqrt(length(x)))
#}

### read kraken report and combine with lineage data ###
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

lineageFile <- "results/oral_diversity/taxa.df"
reports <- list.files("results/oral_diversity/kraken_reports/", full.names=T)

###########################################
############### Wraggling #################
###########################################

# generate lists and bind them to a single df
print("loading reports to memory.")
data_list <- lapply(reports, function(x) readKReport(x, lineageFile))

data_merged <- do.call("rbind", data_list)

grps <- setdiff(names(data_merged), c("perc", "frag", "frag_uniq"))

# use total frag counts across all the files
df_total <- data_merged %>% 
					group_by_at(grps) %>%
					summarise(tot_perc = sum(perc), 
								tot_frag = as.numeric(sum(frag)), 
								tot_frag_uniq = sum(frag_uniq),
								count=n()) %>%
					arrange(desc(tot_perc)) %>%
					ungroup()

# correct column names - issues with taxonkit misassigning at higher taxa levels
	# domian and superkingdoms are under kingdom
	# kingdoms are not present
	# others classifications are correct 
# include only bacteria and virus data

df_corrected <- df_total %>%
					mutate(domain_supking = 
							ifelse(is.na(kingdom), name, as.character(kingdom))) %>%
					dplyr::select(-kingdom) %>%
					filter(domain_supking %in% c("Bacteria", "Viruses"))
# abundance
tot_rank_frag <- df_corrected %>% 
					group_by(rank) %>%
					summarise(tot_rank_frag = sum(tot_frag))
df_abund <- df_corrected %>%
				left_join(tot_rank_frag, by='rank') %>% 
				mutate(rel_abundance = tot_frag/tot_rank_frag)

### species only ###

# bacteria
df_species_bact <- df_abund %>%
					filter(rank %in% c("species")) %>%
					filter(domain_supking == "Bacteria")
		
# virus				
df_species_virus <- df_abund %>%
					filter(rank %in% c("species")) %>%
					filter(domain_supking == "Viruses")

###########################################
############### Analysis ##################
###########################################

### diversity and effective species number (ESN) ###
# bacteria 
bacteria_counts <- df_species_bact$tot_frag
shannon_bact <- diversity(bacteria_counts, index = "shannon")
eff_species_bact <- exp(shannon_bact) # effective number of species

# virus 
virus_counts <- df_species_virus$tot_frag
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

pal <- wes_palette("Darjeeling1", nrow(data_head_bact), type = "continuous")
#pal <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
#			"#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")


#log(tot_frag)
g <- ggplot(data=data_head_bact, aes(x=name, y=tot_frag, fill =name)) +
		facet_grid(vars(phylum), scales="free", space = "free")+
		geom_col(colour='black', show.legend = FALSE)+
		scale_fill_manual(values = pal)+
		theme_classic() +
		ylab("total kmer match") +
		xlab("") +
		coord_flip()+
		theme(strip.text.y = element_text(size = 15, angle=0),
				text = element_text(size=20))

pdf("results/oral_diversity/oralDiv_bacteria.pdf", 15, 20)
#options(scipen=10000)
print(g)
invisible(dev.off())


### virus ###
data_head_virus <- df_species_virus %>% 
					arrange(desc(rel_abundance)) %>%
					head(ceiling(eff_species_virus))

pal <- wes_palette("Darjeeling1", nrow(data_head_virus), type = "continuous")
#pal <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
#			"#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")


#log(tot_frag)
p <- ggplot(data=data_head_virus, aes(x=name, y=tot_frag, fill =name)) +
		facet_grid(vars(phylum), scales="free", space = "free")+
		geom_col(colour='black', show.legend = FALSE)+
		scale_fill_manual(values = pal)+
		theme_classic() +
		ylab("total kmer match") +
		xlab("") +
		coord_flip()+
		theme(strip.text.y = element_text(size = 15, angle=0),
				text = element_text(size=20))

pdf("results/oral_diversity/oralDiv_virus.pdf", 15, 20)
#options(scipen=10000)
print(p)
invisible(dev.off())



system("open -a Skim.app results/oral_diversity/oralDiv_bacteria.pdf")
system("open -a Skim.app results/oral_diversity/oralDiv_virus.pdf")








