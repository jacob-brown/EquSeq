# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-03-18
# Last Modified: 2020-05-29

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

# pavian for shiny app
#pavian::runApp(port=5000)

###########################################
############## Function(s) ################
###########################################

### standard error ###
se <- function(x){
	sd(x)/(sqrt(length(x)))
}

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
data_list <- lapply(reports[1:3], function(x) readKReport(x, lineageFile))

data_merged <- do.call("rbind", data_list)

grps <- setdiff(names(data_merged), c("perc", "frag", "frag_uniq"))

df_summaries <- data_merged %>% 
					group_by_at(grps) %>%
					summarise(avg_perc = mean(perc), 
								se_perc = se(perc), 
								avg_frag = mean(frag), 
								se_frag = se(frag), 
								avg_frag_uniq = mean(frag_uniq),
								se_frag_uniq = se(frag_uniq)) %>%
					arrange(desc(avg_perc)) %>%
					ungroup()



### species only 
data_species <- df_summaries %>%
			filter(rank %in% c("species"))


data_grp <- data_species %>%
				mutate(species = as.character(species),
						kingdom = as.character(kingdom),
						phylum = as.character(phylum)
						) %>%
				mutate(name_grp = 
					ifelse(avg_perc > 10, species, phylum)
					) %>%
				group_by(name_grp, kingdom, phylum) %>%
				summarise(per_grp = sum(avg_perc)) %>%
				filter(per_grp > 0) %>%
				arrange(desc(per_grp)) %>%
				mutate(group_col = 
					ifelse(kingdom != "Eukaryota", kingdom,
						ifelse(phylum == "Chordata", "Animalia",
							ifelse(phylum == "Streptophyta", "Plantae",
								ifelse(phylum == "Ascomycota", "Fungi", NA)))),
					grp=as.factor(1)
					)


# add unclassified/other
data_all <- data_grp %>%
				bind_rows(
					data.frame(name_grp="unclassified/other", 
								per_grp=(100-sum(data_grp$per_grp)))
					) %>%
				mutate(grp=as.factor(1))

# add unclassified
				
data_head <- data_species %>% 
				arrange(desc(avg_perc)) %>%
				group_by(kingdom) %>%
				top_n(10) %>%
				filter(avg_perc > 0.1) %>%
				filter(!(name %in% c("Equus caballus", "Homo sapiens")) &
							phylum != "Streptophyta")


#unclassified <- filter(df_summaries, rank %in% c("unclassified"))
#unq_ranks <- unique(df$rank)
#sum(species$avg_perc) + unclassified$avg_perc


###########################################
############### Plotting ##################
###########################################
pal <- wes_palette("Darjeeling1", nrow(data_head), type = "continuous")
#pal <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
#			"#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")



g <- ggplot(data=data_head, aes(x=name, y=avg_perc, fill =name)) +
		facet_grid(vars(phylum), scales="free", space = "free")+
		geom_col(colour='black', show.legend = FALSE)+
		scale_fill_manual(values = pal)+
		theme_classic() +
		ylab("% reads uniquely identified") +
		xlab("") +
		coord_flip()+
		theme(strip.text.y = element_text(size = 15, angle=0),
				text = element_text(size=20))

pdf("results/oral_diversity/oralDiv.pdf", 9, 7)
print(g)
invisible(dev.off())
system("open -a Skim.app results/oral_diversity/oralDiv.pdf")

#scale_fill_brewer(palette = "Dark2")+


#data_all$name_grp <- factor(data_all$name_grp, levels = rev(data_all$name_grp))
#
#pal <- wes_palette("Darjeeling1", nrow(data_all), type = "continuous")
#pal[1] <- "#B8B8B8"
#
#p <- ggplot(data=data_all, aes(x=grp, y=per_grp, fill = name_grp)) +
#		geom_bar(position="stack", stat="identity", colour="black") +
#		scale_fill_manual(values = pal,
#							name="Taxa")+
#		theme_classic() +
#		theme(axis.title.x=element_blank(),
#				axis.text.x=element_blank(),
#        		axis.ticks.x=element_blank())+
#		ylab("% reads uniquely identified") +
#		scale_y_continuous(expand = c(0,0)) 
# 
#pdf("results/oral_diversity/oralDivFull.pdf", 4, 8)
#print(p)
#invisible(dev.off())
#system("open -a Skim.app results/oral_diversity/oralDivFull.pdf")

