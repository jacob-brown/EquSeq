# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-03-18
# Last Modified: 2020-03-18

# Desc: 


###########################################
################# Modules #################
###########################################

require(ggplot2)

###########################################
############## Function(s) ################
###########################################



###########################################
######### Input(s) and Parameters #########
###########################################

files <- list.files('../results/kraken',  full.names = T)

###########################################
############### Wraggling #################
###########################################

# read in data
df <- read.delim(files[1], sep="\t",  header = F)

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

# filter
species_df <- subset(df, df$rank_code == 'S') # species only
species_df_high <- subset(species_df,  # 0.1% of coverage
							species_df$percentage_cover >= .1)
							
# unique taxon IDs to string
taxonID <- unique(species_df_high$taxon_id)
taxonID_str <- paste(as.character(taxonID), collapse=',')


### get higher level taxa ###
	# bacteria, fungi, etc.
command <- sprintf("efetch -db taxonomy -id %s -format docsum |  
	xtract -pattern DocumentSummary -element  TaxId GenBankDivision", taxonID_str) 

# run command and catch output
out <- system(command, intern = TRUE) 

store <- c() 

# split the strings
for(i in out){
	tmp <- strsplit(i, '\t')[[1]]
	store <- append(store, tmp)
}

# matrix and transpose
t(matrix(store, nrow = 2))

#efetch -db taxonomy -id 3469 -format docsum
#efetch -db taxonomy -id 9796 -format docsum

#Division
###########################################
############### Analysis ##################
###########################################




###########################################
############### Plotting ##################
###########################################



