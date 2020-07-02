# -*- coding: utf-8 -*-
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-06-22
# Last Modified: 2020-06-30

# Desc: 


# first run infer_traits.py

# absent
# heterozygous
# homozygouse
# not present



###########################################
################# Modules #################
###########################################

require(tidyverse)
require(ggplot2)

###########################################
############## Function(s) ################
###########################################



###########################################
######### Input(s) and Parameters #########
###########################################

#df = read.csv("results/gene_to_trait/variants_info.csv")
df = read.csv("results/gene_to_trait/beaglesub.csv")
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

###########################################
############### Wraggling #################
###########################################

# filter(ind %in% c("Ind171", "Ind171.1", "Ind171.2")) %>% # REMOVE ME
#ifelse(str_detect(ind, ".1"), "hetero", 
#								ifelse(str_detect(ind, ".2"), "homo_min", 
#									"homo_maj"))

dft <- df %>%
		dplyr::select(-X) %>%
		pivot_longer(cols = starts_with("Ind"), 
				names_to = "ind", values_to = "GL") %>% 
		mutate(type = ifelse(str_sub(ind, start= -2) == ".1", "hetero", 
								ifelse(str_sub(ind, start= -2) == ".2", "homo_min", 
									"homo_maj")),
				ind = gsub("\\..*","", ind),
				row = row_number()) %>%
		pivot_wider(names_from = type, values_from = GL) %>%
		replace(is.na(.), 0) %>%
		group_by(marker, phen, ref, new, 
					allele1, allele2, ind) %>%
		summarise(homo_maj = sum(homo_maj),
					homo_min = sum(homo_min),
					hetero = sum(hetero)) 

# dft %>% filter(ind == "ind28")
# cat, val
# absent, NA
# homo, GL
# hetero, GL

# return the index of the maximum value
	# draws return both
maxValInd <- function(vec) return(which(max(vec) == vec))

homhet <- c("homo_maj", "homo_min", "hetero")

dft$new <- as.character(dft$new)
dft$allele1 <- as.character(dft$allele1)
dft$allele2 <- as.character(dft$allele2)

tmp <- dft %>%
			mutate(absent = ifelse(new == allele1, FALSE,
								ifelse(new == allele2, FALSE, TRUE))) %>%
			data.frame()


tmp$type <- as.character(NA)
tmp$value <- as.numeric(NA)

# determine highest value category and value
print("determining genotypes.")
for(i in 1:nrow(tmp)){
	
	elem <- tmp[i,]
	type_select <- homhet[maxValInd(elem[homhet])]  # type
	val_select <-  as.vector(elem[type_select][1,]) # GL value
	if(length(type_select) > 1){
		type_select <- "NN"
		val_select <- as.numeric(NA)
	}
	# update
	tmp[i,]$type <- type_select
	tmp[i,]$value <- val_select 
}

# check if major or minor and compare to "new" trait allele
tmp_logic <- tmp %>%
				mutate(new_present = 
					ifelse(absent, FALSE, # if not in pop at all
						ifelse(type == "hetero", TRUE, # hetero must be present
							ifelse(type == "homo_maj",  
								ifelse(new == allele1, TRUE, FALSE), # major
								ifelse(new == allele2, TRUE, FALSE)))), # minor 
						value = ifelse((new_present | absent), value, as.numeric(NA)) # update the value
						) %>%
				mutate(category = 
						ifelse(type == "NN", "ungenotyped", # ensure categorised first
							ifelse(!new_present, NA, # absent 
								ifelse(type == "hetero", "hetero", "homo"))),
						value = ifelse(is.na(category), NA, value)
						) # categories for plot



tmp_logic$ind <- factor(tmp_logic$ind, levels = str_sort(unique(tmp_logic$ind), numeric = T))


p <- ggplot(tmp_logic, mapping = aes(ind, paste(marker,phen), fill = category)) + 
  		geom_tile() +
  		geom_text(aes(label = round(value, 1))) +
    	scale_fill_manual(values = cbbPalette, na.value="white")+
    	theme(axis.text.x = element_text(angle = 90))


pdf("results/gene_to_trait/heat.pdf", 30, 10)
print(p)
dev.off()
system("open -a Skim.app results/gene_to_trait/heat.pdf")




#write.csv(tmp_logic, "results/gene_to_trait/logic.csv")


###########

#tmp_logic %>% filter(ind =="Ind5") %>% data.frame() 
#
#
#
tmp_logic %>%
	filter(ind =="Ind171") %>%
	data.frame()
#	
#tmp_logic %>%
#	filter(marker == "chr3_77739558" & ind =="Ind0")
#
#tmp_logic %>%
#	filter(marker == "chr3_77739558") 



###########################################
############### Analysis ##################
###########################################



###########################################
############### Plotting ##################
###########################################