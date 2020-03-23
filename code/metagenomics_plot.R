# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-03-18
# Last Modified: 2020-03-19

# Desc: 


###########################################
################# Modules #################
###########################################

require(ggplot2)
require(dplyr)
require(tidyr)
#require(metacoder)
require(wesanderson)

###########################################
############## Function(s) ################
###########################################



###########################################
######### Input(s) and Parameters #########
###########################################

#files <- list.files('../results/kraken',  full.names = T)

df <- read.delim('../results/pavian.tsv', sep="\t")
df_all <- read.delim('../results/pavian_all.tsv', sep="\t")

###########################################
############### Wraggling #################
###########################################

# remove unwanted columns
df_fil <- df[,-c(2,3,4, ncol(df))]

lineages <- df[c(1, ncol(df))]

unmapped <- df_all[1,-c(2,3,4, ncol(df))] %>%
				gather("key", "val", -name) %>%
				group_by(name) %>%
				summarise(mean_perc = mean(val), lineage="")


df_new <- gather(df_fil, "key", "val", -name) %>%
			mutate(val=as.numeric(val)) %>% 
			group_by(name) %>%
			summarise(mean_perc = mean(val)) %>%
			arrange(desc(mean_perc)) %>%
			left_join(lineages, by='name') %>%
			dplyr::filter(mean_perc >= .5) 

df_new <- rbind(unmapped, df_new)




df_new$name <- factor(df_new$name, levels = rev(df_new$name))


pal <- wes_palette("Darjeeling1", nrow(df_new), type = "continuous")

pdf("../results/kraken.pdf", 5, 5)

ggplot(data=df_new, aes(x=name, y=mean_perc, fill =name)) +
	geom_col(colour='black', show.legend = FALSE)+
	scale_fill_manual(values = pal)+
	coord_flip() +
	theme_classic() +
	ylab("Percentage Detected") +
	xlab("") +
	ylim(0, 60)

dev.off()















if(F){


					obj <- parse_tax_data(df_new,
					      class_cols = "lineage", # the column that contains taxonomic information
					      class_sep = ">"#, # The character used to separate taxa in the classification
					      #class_regex = "^(.+)__(.+)$", # Regex identifying where the data for each taxon is
					#      class_key = c(tax_rank = "info", # A key describing each regex capture group
					#                    tax_name = "taxon_name")
					)

					# difference in lengths, append 0s
					taxa_n <- length(obj$taxon_names())
					n_zeros <- taxa_n - length(obj$data$tax_data$mean_perc)
					node_col <- c(obj$data$tax_data$mean_perc)

					set.seed(1)

					pdf("../results/kraken.pdf", 5, 5)

					heat_tree(obj, 
					          node_label = obj$taxon_names(),
					          node_size = obj$n_obs(),
					          node_color = node_col, #obj$data$tax_data$mean_perc,
					          node_size_axis_label = "n obs",
					          node_color_axis_label = "Percentage coverage",
					          layout = "davidson-harel", # The primary layout algorithm
					          initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations

					dev.off()





					if(F){
					# storage - create long matrix with appends
					df <- data.frame(matrix(ncol = 6, nrow = 0))

					for(f in files){
						# read in data
						tmp_df <- read.delim(f[1], sep="\t",  header = F)
						df <- rbind(df, tmp_df) # append

					}


					# 1. Percentage of fragments covered by the clade rooted at this taxon
					# 2. Number of fragments covered by the clade rooted at this taxon
					# 3. Number of fragments assigned directly to this taxon
					# 4. A rank code
					# 5. NCBI taxonomic ID number
					# 6. Indented scientific name
					header <- c('percentage_cover', 'number_cover_clade',  'number_cover_taxon',
					 			'rank_code', 'taxon_id', 'name')

					# update the header
					colnames(df) <- header

					df$rank_code <- as.character(df$rank_code)
					df$percentage_cover <- as.numeric(df$percentage_cover)

					# filter
					species_df <- subset(df, df$rank_code == 'S') # species only
					species_df_high <- subset(species_df,  # 0.1% of coverage
												species_df$percentage_cover >= .1)

					species_df_high %>% 
						group_by(name) %>%
						summarise(mean(percentage_cover))
												
					# unique taxon IDs to string
					taxonID <- unique(species_df_high$taxon_id)
					taxonID_str <- paste(as.character(taxonID), collapse=',')


					### get higher level taxa ###
						# bacteria, fungi, etc.
					command <- sprintf("efetch -db taxonomy -id %s -format docsum |  
						xtract -pattern DocumentSummary -element  TaxId GenBankDivision", taxonID_str) 

					#efetch -db taxonomy -id 3469 -format docsum
					#efetch -db taxonomy -id 9796 -format docsum

					#Division

					# run command and catch output
					out <- system(command, intern = TRUE) 

					store <- c() 

					# split the strings
					for(i in out){
						tmp <- strsplit(i, '\t')[[1]]
						store <- append(store, tmp)
					}

					# matrix and transpose
					mat <- t(matrix(store, nrow = 2))
					mat <- cbind(mat[,2], mat[,1])


					#kraken-biom V300044309_L2_B5GHORlfyRAAAAAAA-518.kreport --fmt tsv
					} 







					###########################################
					############### Analysis ##################
					###########################################




					###########################################
					############### Plotting ##################
					###########################################



}